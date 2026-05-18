#!/usr/bin/env python3
"""记录新的平台提交结果"""

import json
import sys
import shutil
from datetime import datetime
from pathlib import Path

SUBMISSIONS_DIR = Path(__file__).parent.parent / "docs" / "submissions"
OUTPUT_DIR = Path(__file__).parent.parent / "output"


def prompt_input(prompt: str, default: str = "") -> str:
    if default:
        value = input(f"{prompt} [{default}]: ").strip()
        return value if value else default
    return input(f"{prompt}: ").strip()


def main():
    print("=" * 60)
    print("记录新的平台提交结果")
    print("=" * 60)
    print()

    # 获取当前 git commit
    try:
        import subprocess
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            text=True,
            cwd=Path(__file__).parent.parent
        ).strip()
    except Exception:
        git_commit = ""

    # 时间
    now = datetime.now()
    default_time = now.strftime("%Y-%m-%d %H:%M:%S")
    time_str = prompt_input("提交时间 (YYYY-MM-DD HH:MM:SS)", default_time)

    # 分数
    score = float(prompt_input("总分数"))
    mol_score = float(prompt_input("mol_score"))
    route_score = float(prompt_input("route_score"))
    sa_score = float(prompt_input("sa_score"))
    validity_score = float(prompt_input("validity_score"))
    binding_score = float(prompt_input("binding_score"))
    route_validity_score = float(prompt_input("route_validity_score"))
    starting_material_score = float(prompt_input("starting_material_availability_score"))
    sample_count = int(prompt_input("sample_count"))
    llm_score = float(prompt_input("llm_score"))

    # 元信息
    exp_round = prompt_input("实验轮次 (如 3)", "")
    hyp_id = prompt_input("假设 ID (如 H001)", "")
    notes = prompt_input("改动说明", "")

    # 解析时间生成文件名
    try:
        dt = datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S")
        filename = dt.strftime("%Y%m%d_%H%M%S")
    except ValueError:
        filename = now.strftime("%Y%m%d_%H%M%S")
        print(f"⚠️  时间格式错误，使用当前时间: {filename}")

    # 备份 result.zip
    result_zip = OUTPUT_DIR / "result.zip"
    backup_zip = OUTPUT_DIR / f"result_{filename}.zip"
    if result_zip.exists():
        shutil.copy2(result_zip, backup_zip)
        print(f"✅ 已备份产物到: {backup_zip}")
    else:
        print(f"⚠️  未找到 result.zip，跳过备份")

    # 生成 JSON
    record = {
        "submission_time": time_str,
        "score": score,
        "breakdown": {
            "mol_score": mol_score,
            "route_score": route_score,
            "sa_score": sa_score,
            "validity_score": validity_score,
            "binding_score": binding_score,
            "route_validity_score": route_validity_score,
            "starting_material_availability_score": starting_material_score,
            "sample_count": sample_count,
            "llm_score": llm_score,
        },
        "experiment_round": int(exp_round) if exp_round else None,
        "hypothesis_id": hyp_id if hyp_id else None,
        "git_commit": git_commit,
        "result_zip_path": str(backup_zip) if result_zip.exists() else "",
        "notes": notes,
    }

    json_path = SUBMISSIONS_DIR / f"{filename}.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=2, ensure_ascii=False)

    print()
    print(f"✅ 提交记录已保存: {json_path}")
    print()

    # 显示最新列表
    print("最新提交列表:")
    subprocess.run([sys.executable, str(Path(__file__).parent / "show_submissions.py")])


if __name__ == "__main__":
    main()
