# 代码演进报告 — Round X (H030)

## 假设
**H030**: LogP-Penalized Composite Scoring for Trivial Route Reduction

## 改动文件

### `tools/pipeline.py`

#### 改动1：新增 `compute_logp_score()` 函数（第46-67行）
```python
def compute_logp_score(logp: float) -> float:
    """Compute LogP reasonableness score (H030)."""
    if logp <= 2.0:
        return 1.0
    if logp >= 5.0:
        return 0.0
    return max(0.0, 1.0 - (logp - 2.0) / 3.0)
```

**设计理由**: LogP≤2 为理想药物样分子（满分1.0），LogP≥5 为 Lipinski 违规（0分），中间线性衰减。这样 LogP=3.5 的多环芳烃得分 0.5，LogP=2.5 的铰链结合分子得分 0.83。

#### 改动2：候选分子记录新增 `logp` 字段（第256行）
```python
"logp": mol.get("logp", 0.0),
```

#### 改动3：H012 复合评分公式修改（第277-296行）
```
旧: 0.80×BE_norm + 0.20×route_quality
新: 0.75×BE_norm + 0.15×route_quality + 0.10×logp_score
```

**权重设计**: BE 从 80% 降至 75%，路线质量从 20% 降至 15%，新增 LogP 10%。预期的净效应：
- LogP=2.0 分子：logp_score=1.0 → +0.10 加成
- LogP=3.5 分子：logp_score=0.5 → +0.05 加成
- LogP=5.0 分子：logp_score=0.0 → +0.00 加成

对高 BE 但高 LogP（trivial 倾向）的分子，logp_score 低将抵消其 BE 优势。

#### 改动4：日志输出显示 logP（第305-308行）

## 兼容性
- 保留所有现有接口不变
- `compute_logp_score` 是 pipeline.py 内部函数
- 不修改 `src/scorer.py`（`compute_total_score` 仍用于事后评估）
- 不修改进化过程（`generate_with_docking_guidance` 中排序保持纯 BE）

## 编译验证
```
$ python3 -m py_compile tools/pipeline.py
OK ✅
```
