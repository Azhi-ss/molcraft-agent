#!/usr/bin/env python3
"""解析 result.log 并转换成结构化 JSON 格式。"""

import json
import re
from pathlib import Path
from typing import Any


def parse_log(log_path: Path) -> dict[str, Any]:
    """解析日志文件并返回结构化字典。"""
    content = log_path.read_text(encoding="utf-8")
    lines = content.split("\n")

    result: dict[str, Any] = {
        "metadata": {},
        "stages": [],
        "docking_progress": {},
        "top_molecules": [],
        "experiment_results": {},
        "hypothesis_validation": {},
        "iteration_report": {}
    }

    # 1. 解析元数据
    start_time_match = re.search(r"\[(.*?)\] MolCraft Agent 启动", content)
    if start_time_match:
        result["metadata"]["start_time"] = start_time_match.group(1)

    end_time_match = re.search(r"\[([\d\-T:.]+)\] MolCraft Agent 执行结束", content)
    if end_time_match:
        result["metadata"]["end_time"] = end_time_match.group(1)

    program_match = re.search(r"指令书: (.*?)\n", content)
    if program_match:
        result["metadata"]["program_file"] = program_match.group(1)

    iter_log_match = re.search(r"迭代记录: (.*?)\n", content)
    if iter_log_match:
        result["metadata"]["iteration_log"] = iter_log_match.group(1)

    config_match = re.search(r"配置: 最大 (\d+) 次迭代 \| 最长 (\d+) 分钟", content)
    if config_match:
        result["metadata"]["max_iterations"] = int(config_match.group(1))
        result["metadata"]["max_minutes"] = int(config_match.group(2))

    # 2. 解析实验结果
    avg_match = re.search(r"Top 10 average binding energy: \*\*(.*?) kcal/mol\*\*", content)
    if avg_match:
        result["experiment_results"]["top10_avg_binding_energy"] = float(avg_match.group(1))

    best_match = re.search(r"Best binding energy: \*\*(.*?) kcal/mol\*\*", content)
    if best_match:
        result["experiment_results"]["best_binding_energy"] = float(best_match.group(1))

    trivial_match = re.search(r"Trivial route ratio: \*\*(.*?)\*\*", content)
    if trivial_match:
        result["experiment_results"]["trivial_route_ratio"] = trivial_match.group(1)

    success_match = re.search(r"Docking success rate: (.*?) ", content)
    if success_match:
        result["experiment_results"]["docking_success_rate"] = success_match.group(1)

    # 历史基线对比
    baseline_avg_match = re.search(r"Top 10 avg: (.*?) kcal/mol", content)
    if baseline_avg_match:
        result["experiment_results"]["baseline_top10_avg"] = float(baseline_avg_match.group(1))

    baseline_best_match = re.search(r"Best: (.*?) kcal/mol", content)
    if baseline_best_match:
        result["experiment_results"]["baseline_best_binding_energy"] = float(baseline_best_match.group(1))

    # 3. 解析分子对接进度
    docking_matches = re.findall(r"\[对接\] (\d+)/(\d+):", content)
    if docking_matches:
        last = docking_matches[-1]
        result["docking_progress"]["total"] = int(last[1])
        result["docking_progress"]["completed"] = int(last[0])
        result["docking_progress"]["success_rate"] = f"{int(last[0])/int(last[1])*100:.1f}%"

    # 4. 从 result.csv 输出先建立完整的 10 个分子列表
    csv_pattern = r"mol_smiles,route\n(.*?)\n\n"
    csv_match = re.search(csv_pattern, content, re.DOTALL)
    if csv_match:
        csv_lines = csv_match.group(1).strip().split("\n")
        result["top_molecules"] = []
        for line in csv_lines:
            if ">>" in line:
                parts = line.split(",", 1)
                if len(parts) == 2:
                    smiles, route = parts
                    result["top_molecules"].append({
                        "smiles": smiles,
                        "synthesis_route": route
                    })

    # 5. 从合成路线输出补充 steps 和 trivial 信息
    route_pattern = r"(\S+)\s+steps=(\d+) trivial=(True|False) route=(\S+)"
    route_matches = re.findall(route_pattern, content)
    for smiles, steps, trivial, route in route_matches:
        for mol in result["top_molecules"]:
            if mol["smiles"] == smiles:
                mol["synthesis_steps"] = int(steps)
                mol["trivial_route"] = trivial == "True"
                break

    # 6. 从 SA score 对比输出补充 SA 和 rings 信息
    sa_matches = re.findall(r"SA=([\d.]+) rings=(\d+) pass_new=(\w+) pass_old=(\w+)\s+(\S+)", content)
    for sa, rings, pass_new, pass_old, smiles in sa_matches:
        for mol in result["top_molecules"]:
            if mol["smiles"] == smiles:
                mol["sa_score"] = float(sa)
                mol["rings"] = int(rings)
                mol["pass_new_filter"] = pass_new == "True"
                mol["pass_old_filter"] = pass_old == "True"
                break

    # 6.1 从 docs/experiment_round_3.md 提取完整的 top-10 分子数据（包括结合能）
    report_path = Path("docs/experiment_round_3.md")
    if report_path.exists():
        report_content = report_path.read_text(encoding="utf-8")
        # 提取表格中的分子数据 - 处理 "0 (trivial)" 格式的步数
        table_pattern = r"\|\s+(\d+)\s+\|\s+(\S+)\s+\|\s+([-\d.]+)\s+\|\s+([\d.]+)\s+\|\s+(\d+)\s+\|\s+([^|]+)\s+\|\s+([^|]+)\s+\|"
        table_matches = re.findall(table_pattern, report_content)
        for idx, smiles, be, sa, rings, steps_str, quality in table_matches:
            steps_str = steps_str.strip()
            is_trivial = "trivial" in steps_str.lower() or steps_str == "0"
            steps = 0 if is_trivial else int(steps_str) if steps_str.isdigit() else None

            for mol in result["top_molecules"]:
                if mol["smiles"] == smiles:
                    mol["binding_energy"] = float(be)
                    mol["sa_score"] = float(sa)
                    mol["rings"] = int(rings)
                    mol["synthesis_steps"] = steps
                    mol["route_quality"] = quality.strip()
                    mol["trivial_route"] = is_trivial
                    break

    # 7. 解析假设验证
    hypo_match = re.search(r"假设 (H\d+) 验证成功", content)
    if hypo_match:
        result["hypothesis_validation"]["hypothesis_id"] = hypo_match.group(1)
        result["hypothesis_validation"]["success"] = True

    summary_match = re.search(r"SA score 阈值收紧（(\d+\.\d+)→(\d+\.\d+)）\+ 环数限制（≤(\d+)）", content)
    if summary_match:
        result["hypothesis_validation"]["changes"] = {
            "sa_threshold_old": float(summary_match.group(1)),
            "sa_threshold_new": float(summary_match.group(2)),
            "max_rings": int(summary_match.group(3))
        }

    # 8. 解析迭代报告
    round_match = re.search(r"第 (\d+) 轮迭代已记录", content)
    if round_match:
        result["iteration_report"]["round"] = int(round_match.group(1))
        result["iteration_report"]["status"] = "completed"

    report_match = re.search(r"实验报告已写入 `(.*?)`", content)
    if report_match:
        result["iteration_report"]["report_file"] = report_match.group(1)

    # 8. 解析完整的阶段内容
    stage_pattern = r"## (.*?)\n(.*?)(?=\n## |\n={50,}|\n\[2026)"
    stage_matches = re.findall(stage_pattern, content, re.DOTALL)
    for stage_name, stage_content in stage_matches:
        stage_name = stage_name.strip()

        # 提取系统命令输出
        system_pattern = r"<system>(.*?)</system>"
        system_commands = re.findall(system_pattern, stage_content, re.DOTALL)

        # 清理内容（移除系统命令标记）
        clean_content = re.sub(r"<system>.*?</system>", "", stage_content, flags=re.DOTALL).strip()

        result["stages"].append({
            "name": stage_name,
            "content": clean_content,
            "system_commands_count": len(system_commands)
        })

    # 9. 解析下一步建议
    next_steps_match = re.search(r"可以考虑以下假设之一：(.*?)您希望继续还是就此完成？", content, re.DOTALL)
    if next_steps_match:
        next_steps_text = next_steps_match.group(1).strip()
        suggestions = re.findall(r"\d+\. \*\*(H\d+)\*\* — (.*?)(?=\n\d+\. |$)", next_steps_text, re.DOTALL)
        result["next_suggestions"] = [
            {"hypothesis_id": h, "description": d.strip()}
            for h, d in suggestions
        ]

    return result


def main() -> None:
    log_path = Path("output/result.log")
    if not log_path.exists():
        print(f"错误: 找不到 {log_path}")
        return

    print(f"正在解析: {log_path}")
    result = parse_log(log_path)

    # 输出 JSON
    output_path = Path("output/result.json")
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    print(f"已生成: {output_path}")
    print(f"文件大小: {output_path.stat().st_size / 1024:.1f} KB")

    # 打印摘要
    print("\n=== 解析摘要 ===")
    meta = result["metadata"]
    print(f"运行时间: {meta.get('start_time', 'N/A')} → {meta.get('end_time', 'N/A')}")
    print(f"执行阶段: {len(result['stages'])} 个")
    print(f"分子数量: {len(result['top_molecules'])}")

    res = result["experiment_results"]
    print(f"\n=== 实验结果 ===")
    print(f"Top-10 平均结合能: {res.get('top10_avg_binding_energy', 'N/A')} kcal/mol")
    print(f"最优结合能: {res.get('best_binding_energy', 'N/A')} kcal/mol")
    print(f"Trivial 路线比例: {res.get('trivial_route_ratio', 'N/A')}")
    print(f"对接成功率: {res.get('docking_success_rate', 'N/A')}")

    hypo = result["hypothesis_validation"]
    print(f"\n=== 假设验证 ===")
    print(f"假设ID: {hypo.get('hypothesis_id', 'N/A')}")
    print(f"验证结果: {'✅ 成功' if hypo.get('success') else '❌ 失败'}")

    if "changes" in hypo:
        print(f"参数变更:")
        print(f"  SA 阈值: {hypo['changes']['sa_threshold_old']} → {hypo['changes']['sa_threshold_new']}")
        print(f"  最大环数: ≤{hypo['changes']['max_rings']}")


if __name__ == "__main__":
    main()
