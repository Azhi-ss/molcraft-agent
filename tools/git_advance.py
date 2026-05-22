#!/usr/bin/env python3
"""Git 实验管理脚本 —— 参考 autoresearch 的 advance/revert 模式。

每次实验迭代后：
- 成功（结合能改善）→ git commit 保留
- 失败（结合能退化）→ git reset 回退

用法:
    python tools/git_advance.py --round 1 --best-be -7.92 --status keep
    python tools/git_advance.py --round 2 --best-be -7.50 --status discard
"""
import argparse
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str]) -> str:
    """运行 shell 命令并返回 stdout。"""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {' '.join(cmd)}\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result.stdout.strip()


def main():
    parser = argparse.ArgumentParser(description="Git 实验管理")
    parser.add_argument("--round", type=int, required=True, help="实验轮次")
    parser.add_argument("--best-be", type=float, required=True, help="本轮最佳结合能")
    parser.add_argument(
        "--status",
        choices=["keep", "discard", "crash"],
        required=True,
        help="keep=保留 commit, discard=回退到上一轮, crash=回退并记录失败",
    )
    parser.add_argument(
        "--files",
        nargs="*",
        default=None,
        help="discard 时精确还原的文件列表（可选，不指定则还原全部）",
    )
    args = parser.parse_args()

    # 检查是否在 git 仓库中
    try:
        run(["git", "rev-parse", "--git-dir"])
    except SystemExit:
        print("错误：当前目录不是 git 仓库。请先初始化 git。")
        sys.exit(1)

    if args.status == "keep":
        # 保留本轮实验结果
        msg = f"round-{args.round}: best_BE={args.best_be:.3f} kcal/mol"
        run(["git", "add", "-A"])
        run(["git", "commit", "-m", msg, "--allow-empty"])
        print(f"✅ 已提交: {msg}")
        print(f"   commit: {run(['git', 'rev-parse', '--short', 'HEAD'])}")

    elif args.status in ("discard", "crash"):
        # 回退本轮代码改动
        if args.files:
            # 精确还原指定文件（从 backup commit 恢复）
            for f in args.files:
                try:
                    run(["git", "checkout", "HEAD~1", "--", f])
                    print(f"   已还原: {f}")
                except SystemExit:
                    print(f"   ⚠️ 还原失败: {f}（文件可能在 HEAD~1 不存在）")
            # 还原后提交
            run(["git", "add", "-A"])
            run(["git", "commit", "-m", f"round-{args.round}: discarded (reverted {len(args.files)} files)"])
            print(f"🔄 已精确还原 round-{args.round} 的 {len(args.files)} 个文件")
        else:
            # 全量回退：stash + soft reset
            run(["git", "stash", "push", "-u", "-m", f"round-{args.round}-discard"])
            run(["git", "reset", "--soft", "HEAD~1"])
            status = "丢弃" if args.status == "discard" else "崩溃回退"
            print(f"🔄 已{status} round-{args.round} 的全部更改，回到上一轮状态")


if __name__ == "__main__":
    main()