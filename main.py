#!/usr/bin/env python3
"""MolCraft Agent 入口 —— 使用 Kimi Agent SDK 启动自主药物研发智能体。

本脚本遵循 autoresearch 架构模式：
- program.md：人类编辑的 Agent 指令书（唯一人机接口）
- docs/iteration_log.jsonl：Agent 主动报告的迭代记录
- main.py：只负责读取 program.md 并启动 Agent，自身不携带业务逻辑

运行方式:
    cd molcraft-agent
    source .venv/bin/activate
    python main.py                    # 默认 1 次迭代，最长 90 分钟
    python main.py --iterations 3     # 迭代 3 次
    python main.py --max-minutes 60   # 最长 60 分钟
    python main.py --max-steps 1000   # 每轮最多 1000 步
"""
import argparse
import asyncio
import json
import uuid
import zipfile
import os
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Any

from kaos.path import KaosPath

from src.event_logger import EventLogger
from src.event_schema import MetricsEvent

# 必须在导入 kimi_agent_sdk 之前加载环境变量！
from dotenv import load_dotenv
from pydantic import SecretStr

from kimi_cli.config import Config, LLMModel, LLMProvider

load_dotenv()

# 在导入 SDK 前注入 DeepSeek reasoning_key 支持
import kosong.contrib.chat_provider.openai_legacy as _openai_provider
_original_init = _openai_provider.OpenAILegacy.__init__

def _patched_openai_init(self, *, model, **kwargs):
    # DeepSeek V4 使用 reasoning_content 作为推理内容字段名
    if "deepseek" in model.lower():
        kwargs["reasoning_key"] = "reasoning_content"
    _original_init(self, model=model, **kwargs)

_openai_provider.OpenAILegacy.__init__ = _patched_openai_init

from kimi_agent_sdk import prompt


PROJECT_ROOT = Path(__file__).parent

AGENT_YAML = Path(__file__).parent / "agent.yaml"
PROGRAM_MD = Path(__file__).parent / "program.md"
OUTPUT_DIR = Path(__file__).parent / "output"
ITERATION_LOG = Path(__file__).parent / "docs" / "iteration_log.jsonl"
CUSTOM_LLM_PROVIDER_KEY = "env-openai-compatible"
CUSTOM_LLM_MODEL_KEY = "env-llm"
PIPELINE_CSV_PATH_ENV = "MOLCRAFT_CSV_PATH"


def _parse_int_env(name: str, default: int) -> int:
    value = os.getenv(name)
    if not value:
        return default
    try:
        return int(value)
    except ValueError:
        return default


def build_agent_llm_config() -> tuple[Config | None, str | None]:
    """Build an explicit SDK config from LLM_* env vars when provided.

    kimi_agent_sdk treats ``model=`` as a config model key, not a raw provider
    model name. Supplying an in-memory config keeps LLM_MODEL tied to the
    requested OpenAI-compatible provider instead of falling back to Kimi config.
    """
    base_url = os.getenv("LLM_BASE_URL")
    api_key = os.getenv("LLM_API_KEY")
    if not base_url or not api_key:
        return None, None

    model_name = os.getenv("LLM_MODEL", "deepseek-chat")
    provider_type = os.getenv("LLM_PROVIDER_TYPE", "openai_legacy")
    max_context_size = _parse_int_env("LLM_MAX_CONTEXT_SIZE", 100_000)

    config = Config(
        default_model=CUSTOM_LLM_MODEL_KEY,
        models={
            CUSTOM_LLM_MODEL_KEY: LLMModel(
                provider=CUSTOM_LLM_PROVIDER_KEY,
                model=model_name,
                max_context_size=max_context_size,
            )
        },
        providers={
            CUSTOM_LLM_PROVIDER_KEY: LLMProvider(
                type=provider_type,  # type: ignore[arg-type]
                base_url=base_url,
                api_key=SecretStr(api_key),
            )
        },
    )
    return config, CUSTOM_LLM_MODEL_KEY


def load_program() -> str:
    """读取 program.md，这是人类与 Agent 交互的唯一接口。"""
    if not PROGRAM_MD.exists():
        raise FileNotFoundError(
            f"{PROGRAM_MD} 不存在。这是 Agent 的指令书，必须创建后才能运行。"
        )
    return PROGRAM_MD.read_text(encoding="utf-8")


def zip_results(output_dir: Path) -> None:
    """打包结果文件供平台提交。原子操作：先写临时文件，成功后替换。"""
    result_csv = output_dir / "result.csv"
    result_log = output_dir / "result.log"

    # 使用唯一临时文件名，避免并发冲突
    fd, temp_zip_str = tempfile.mkstemp(prefix="result.", suffix=".zip", dir=output_dir)
    os.close(fd)
    temp_zip = Path(temp_zip_str)

    try:
        with zipfile.ZipFile(temp_zip, "w", zipfile.ZIP_DEFLATED) as zf:
            if result_csv.exists():
                zf.write(result_csv, arcname="result.csv")
            if result_log.exists():
                zf.write(result_log, arcname="result.log")

        # 验证创建的 zip 文件有效
        with zipfile.ZipFile(temp_zip, "r") as zf:
            files = zf.namelist()
            if "result.csv" not in files or "result.log" not in files:
                raise RuntimeError(f"Zip missing required files: {files}")

        # 原子替换
        os.replace(str(temp_zip), str(output_dir / "result.zip"))
        print(f"\n[打包完成] {output_dir / 'result.zip'}")
    except Exception as e:
        if temp_zip.exists():
            temp_zip.unlink(missing_ok=True)
        raise RuntimeError(f"打包失败: {e}") from e


def count_iterations() -> int:
    """统计 Agent 已报告的迭代次数。

    Agent 每完成一轮「文献→诊断→改代码→实验验证」闭环后，
    会主动调用 report_iteration 工具，写入一条记录到 docs/iteration_log.jsonl。
    """
    if not ITERATION_LOG.exists():
        return 0
    count = 0
    try:
        with open(ITERATION_LOG, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    count += 1
    except Exception:
        pass
    return count


async def monitor(stop_event: asyncio.Event, max_iterations: int, max_seconds: float, initial_iteration_count: int):
    """独立协程：每 30 秒检查迭代次数和时间限制。"""
    start_time = time.time()
    limit_printed = False

    while not stop_event.is_set():
        await asyncio.sleep(30)

        elapsed = time.time() - start_time

        # 时间上限
        if elapsed > max_seconds:
            print(
                f"\n\n[{datetime.now().isoformat()}] 达到最大运行时间 "
                f"({int(max_seconds // 60)} 分钟)，正在退出...",
                flush=True,
            )
            stop_event.set()
            return

        # 迭代次数上限（基于本轮新增的迭代，而非累积）
        if max_iterations > 0 and not limit_printed:
            current = count_iterations()
            per_run_count = current - initial_iteration_count
            print(
                f"\n[{datetime.now().isoformat()}] 迭代进度: {per_run_count}/{max_iterations} (初始: {initial_iteration_count})",
                flush=True,
            )
            if per_run_count >= max_iterations:
                print(
                    f"\n[{datetime.now().isoformat()}] 已达到最大迭代次数 "
                    f"({max_iterations} 次)，60 秒后退出...",
                    flush=True,
                )
                limit_printed = True
                await asyncio.sleep(60)
                stop_event.set()
                return


async def main() -> None:
    parser = argparse.ArgumentParser(description="MolCraft Agent 启动器")
    parser.add_argument(
        "--iterations",
        type=int,
        default=1,
        help="最大迭代次数。默认 1 次。设为 0 表示不限制。",
    )
    parser.add_argument(
        "--max-minutes",
        type=int,
        default=90,
        help="最大运行时间（分钟）。默认 90 分钟。",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=1000,
        help="每轮最大步数。默认 1000 步。",
    )
    parser.add_argument(
        "--thinking",
        action="store_true",
        default=True,
        help="启用 DeepSeek reasoning 模式（默认开启）。",
    )
    parser.add_argument(
        "--no-thinking",
        action="store_false",
        dest="thinking",
        help="禁用 DeepSeek reasoning 模式。",
    )
    args = parser.parse_args()

    # 记录本次运行开始时间，用于检测新生成的产物
    run_start_time = time.time()

    # 生成唯一 run_id，用于 evomap 多 session 区分
    run_id = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"
    os.environ["MOLCRAFT_RUN_ID"] = run_id

    program = load_program()
    llm_config, llm_model_key = build_agent_llm_config()

    # 确保输出目录存在
    OUTPUT_DIR.mkdir(exist_ok=True)
    log_path = OUTPUT_DIR / "result.log"

    # 重定向 stdout 到终端 + 文件
    tee = EventLogger(log_path)
    original_stdout = sys.stdout
    sys.stdout = tee

    pipeline_csv_path = tee.tmp_path.with_suffix(".pipeline.tmp.csv")
    pipeline_csv_path.unlink(missing_ok=True)
    os.environ[PIPELINE_CSV_PATH_ENV] = str(pipeline_csv_path.resolve())

    stop_event = asyncio.Event()
    max_seconds = args.max_minutes * 60

    # 记录本次运行前的初始迭代次数，用于计算本轮新增迭代
    initial_iteration_count = count_iterations()

    # 启动独立监控协程
    monitor_task = asyncio.create_task(
        monitor(stop_event, args.iterations, max_seconds, initial_iteration_count)
    )

    try:
        run_status = "success"
        print(f"[{datetime.now().isoformat()}] MolCraft Agent 启动")
        print("=" * 60)
        print("模式: autoresearch（LLM 自主迭代实验）")
        print(f"指令书: {PROGRAM_MD}")
        print(f"迭代记录: {ITERATION_LOG}")
        print(f"日志文件: {log_path}")
        print(f"配置: 最大 {args.iterations} 次迭代 | 最长 {args.max_minutes} 分钟")
        if llm_config and llm_model_key:
            provider = llm_config.providers[CUSTOM_LLM_PROVIDER_KEY]
            model = llm_config.models[llm_model_key]
            print(f"[配置] 使用自定义 LLM: {model.model}")
            print(f"[配置] Provider: {provider.type}")
            print(f"[配置] API Endpoint: {provider.base_url}")
        print("(yolo 模式下自动批准所有操作)")
        print()

        # 记录启动事件
        tee.log_start(round=count_iterations() + 1, max_minutes=args.max_minutes,
                      max_iterations=args.iterations, program_file=str(PROGRAM_MD))

        async for msg in prompt(
            program,
            config=llm_config,
            agent_file=AGENT_YAML,
            skills_dir=KaosPath(PROJECT_ROOT / ".kimi" / "skills"),
            yolo=True,
            thinking=args.thinking,
            max_steps_per_turn=args.max_steps,
            model=llm_model_key,
        ):
            if stop_event.is_set():
                break

            text = msg.extract_text()
            if text:
                print(text, end="", flush=True)

        print()
        print()
        print("=" * 60)
        print(f"[{datetime.now().isoformat()}] MolCraft Agent 执行结束")
        print("请检查:")
        print("  - output/result.csv（最终候选分子与合成路线）")
        print(f"  - {log_path}（详细结果日志）")
        print(f"  - {ITERATION_LOG}（迭代记录）")
        print("=" * 60)
    except Exception:
        run_status = "failure"
        raise
    finally:
        stop_event.set()
        if not monitor_task.done():
            monitor_task.cancel()
            try:
                await monitor_task
            except asyncio.CancelledError:
                pass
        sys.stdout = original_stdout

        # Success validation and stale artifact detection
        result_csv = OUTPUT_DIR / "result.csv"
        if not (pipeline_csv_path.exists() and pipeline_csv_path.stat().st_mtime > run_start_time):
            run_status = "failure"
            print(f"\n检测失败：未发现本轮运行产生的临时 result.csv", flush=True)
        else:
            # 计算 metrics（在 commit 之前，因为日志还写在 tmp 文件）
            try:
                import csv
                molecules = []
                with open(pipeline_csv_path, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        molecules.append(row)

                trivial_count = sum(
                    1 for m in molecules
                    if m.get('route', '') == m.get('mol_smiles', '') + '>>' + m.get('mol_smiles', '')
                )
                non_trivial_count = len(molecules) - trivial_count
                tee.log_metrics(MetricsEvent(
                    molecule_count=len(molecules),
                    non_trivial_count=non_trivial_count,
                    trivial_count=trivial_count,
                    avg_binding_energy=None,
                    min_binding_energy=None,
                    avg_syn_steps=None,
                    docking_success_rate=0.0,
                ))
            except Exception:
                pass

        # 记录结束事件
        tee.log_end(status=run_status)

        # 只有成功运行且有新产物时才提交日志并打包结果
        if run_status == "success":
            # CSV 内容验证：检查 result.csv 格式和内容
            try:
                import csv
                with open(pipeline_csv_path, "r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    rows = list(reader)
                    if len(rows) == 0:
                        raise RuntimeError("result.csv 为空，没有分子数据")
                    if "mol_smiles" not in reader.fieldnames or "route" not in reader.fieldnames:
                        raise RuntimeError(f"result.csv 缺少必要列: {reader.fieldnames}")
                    for row in rows:
                        if not row.get("mol_smiles", "").strip():
                            raise RuntimeError("result.csv 存在空分子行")
                        if not row.get("route", "").strip():
                            raise RuntimeError("result.csv 存在空合成路线行")
                print(f"\nresult.csv 验证通过，共 {len(rows)} 个分子", flush=True)
            except Exception as e:
                run_status = "failure"
                print(f"\nresult.csv 验证失败: {e}", flush=True)

        if run_status == "success":
            # JSONL 格式验证：在 commit 之前验证 tmp 文件
            # 使用 parse_constant 严格禁止 NaN、Infinity、-Infinity
            def _reject_nan_constant(val: str) -> None:
                raise ValueError(f"无效 JSON 数值常量: {val}")

            try:
                with open(tee.tmp_path, "r", encoding="utf-8") as f:
                    for line_num, line in enumerate(f, 1):
                        if line.strip():
                            json.loads(line, parse_constant=_reject_nan_constant)
            except (json.JSONDecodeError, ValueError) as e:
                run_status = "failure"
                print(f"\n日志格式验证失败，行 {line_num}: {e}", flush=True)
            except Exception as e:
                run_status = "failure"
                print(f"\n日志格式验证错误: {e}", flush=True)

            if run_status == "success":
                print(f"\nJSONL 格式验证通过", flush=True)
                os.replace(pipeline_csv_path, result_csv)
                tee.commit()
                zip_results(OUTPUT_DIR)
            # 失败时保留旧的 result.zip，不破坏上一次有效提交

        if run_status != "success" and pipeline_csv_path.exists():
            pipeline_csv_path.unlink(missing_ok=True)
        tee.close()


if __name__ == "__main__":
    asyncio.run(main())
