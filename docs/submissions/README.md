# 平台提交记录

记录每次提交到 MolCraft 平台后的分数和产物。

## 目录结构

```
docs/submissions/
├── README.md
└── YYYYMMDD_HHMMSS.json    # 每次提交的元数据
output/
├── result.zip              # 最新提交产物
└── result_YYYYMMDD_HHMMSS.zip  # 历史提交产物备份
```

## 字段说明

每个 JSON 文件包含：

| 字段 | 说明 |
|------|------|
| `submission_time` | 提交时间 |
| `score` | 平台总分 |
| `breakdown` | 各细分项分数 |
| `experiment_round` | 对应的实验轮次 |
| `hypothesis_id` | 验证的假设 ID |
| `git_commit` | 代码提交 hash |
| `result_zip_path` | 产物备份路径 |
| `notes` | 改动说明 |

## 快速查看

```bash
# 查看所有提交分数
python tools/show_submissions.py
```

## 提交流程

1. 运行 `python main.py` 或 `python tools/pipeline.py`
2. 上传 `output/result.zip` 到平台
3. 得到分数后运行记录脚本
4. 提交代码并记录 git commit hash
