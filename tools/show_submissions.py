#!/usr/bin/env python3
"""快速查看所有平台提交记录"""

import json
import glob
from pathlib import Path
from datetime import datetime

SUBMISSIONS_DIR = Path(__file__).parent.parent / "docs" / "submissions"


def main():
    files = sorted(glob.glob(str(SUBMISSIONS_DIR / "*.json")))

    if not files:
        print("暂无提交记录")
        return

    print(f"{'=' * 100}")
    print(f"{'提交时间':<20} {'总分':>10} {'结合能':>10} {'合成':>10} {'样本':>8} {'假设':>8}")
    print(f"{'=' * 100}")

    total_score = 0.0
    for f in files:
        with open(f) as fp:
            data = json.load(fp)

        time_str = datetime.fromisoformat(data["submission_time"]).strftime("%Y-%m-%d %H:%M")
        score = f"{data['score']:.4f}"
        binding = f"{data['breakdown']['binding_score']:.4f}"
        route = f"{data['breakdown']['route_score']:.4f}"
        samples = f"{data['breakdown']['sample_count']}"
        hyp = data.get("hypothesis_id", "-")

        print(f"{time_str:<20} {score:>10} {binding:>10} {route:>10} {samples:>8} {hyp:>8}")
        total_score += data["score"]

    print(f"{'=' * 100}")
    print(f"共 {len(files)} 次提交，平均总分: {total_score / len(files):.4f}")
    print()


if __name__ == "__main__":
    main()
