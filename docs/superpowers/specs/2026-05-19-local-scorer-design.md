# Local Scoring System Design

## Problem

比赛评分系统的归一化函数未知，导致本地分数与比赛分数偏差大。我们只知道评分公式框架，但 binding_score、sa_score 等子项如何从原始值归一化到 0-1 是黑盒。

## Goal

构建本地评分系统，通过多次提交校准归一化函数，使本地预测分数与比赛实际分数误差 < 5%。

## Architecture

### Components

1. **`src/scorer.py`** — 本地评分引擎
   - 输入：`result.csv`（mol_smiles, route）
   - 输出：完整评分报告（各子项分数 + 总分）
   - 包含参数化归一化函数，校准后更新参数

2. **`src/calibrator.py`** — 校准器
   - 输入：本地 Vina 分数 + 比赛 binding_score 数据对
   - 输出：拟合的归一化函数参数
   - 持久化校准结果到 `data/calibration.json`

3. **`scripts/score_local.py`** — CLI 入口
   - `python scripts/score_local.py output/result.csv` — 评分
   - `python scripts/score_local.py --calibrate` — 用已有数据校准

### Scoring Formula (Known Framework)

```
total_score = 0.7 * mol_score + 0.3 * route_score

mol_score = 0.8 * binding_score
          + 0.1 * validity_score
          + 0.1 * sa_score

route_score = 0.55 * route_validity_score
            + 0.30 * starting_material_availability_score
            + 0.05 * step_penalty_score
            + 0.05 * convergence_score
            + 0.05 * balance_score
```

### Normalization Functions (To Calibrate)

#### binding_score

本地 Vina 输出原始结合能（kcal/mol，负值），需归一化到 0-1。

假设候选函数形式（按优先级）：
1. **Min-Max**: `score = (vina - min) / (max - min)` — 最简单
2. **Sigmoid**: `score = 1 / (1 + exp(-k * (vina - midpoint)))`
3. **Clipped Linear**: `score = clamp((vina - threshold) / range, 0, 1)`

校准数据：本地 Vina 分数 → 比赛 binding_score

#### sa_score

已知：SAScore > 4 → 0，< 4 时越低越好。具体归一化待定。

候选：
1. **Step**: `score = 0 if sa > 4 else (4 - sa) / 4`
2. **Clipped**: `score = max(0, (4 - sa) / 4)`

#### route_validity_score

已知逻辑：每条路线通过/不通过，取通过比例。零化条件明确。

#### starting_material_availability_score

已知：DB hit → 1.0，否则用 SAScore fallback。具体 fallback 逻辑待定。

#### balance_score, step_penalty_score, convergence_score

这些子项权重仅 0.05，对总分影响小。先实现基本逻辑，后续精调。

### Calibration Workflow

```
Phase 1: 构建 + 基线校准
  ├─ 实现 scorer.py 全部子项
  ├─ 用默认参数跑一次本地评分
  ├─ 提交 result.csv 获取比赛分数
  ├─ 对比本地 vs 比赛，识别偏差
  └─ 用偏差数据初步校准归一化参数

Phase 2: 定向校准 binding_score
  ├─ 设计 3 组校准分子（高/中/低 Vina）
  ├─ 分别提交，收集 binding_score
  ├─ 拟合 vina → binding_score 映射
  └─ 更新 calibration.json

Phase 3: 验证 + 生产
  ├─ 用校准后的评分系统筛选分子
  ├─ 本地预测分数 vs 实际提交分数 < 5% 误差
  └─ 用精准评分指导分子优化
```

### Calibration Data Format

```json
{
  "version": 1,
  "calibrated_at": "2026-05-19T...",
  "binding_score": {
    "function": "clipped_linear",
    "params": {
      "threshold": -5.0,
      "range": 10.0
    },
    "calibration_points": [
      {"vina_raw": -9.591, "actual_score": 0.1766},
      {"vina_raw": -8.0, "actual_score": null}
    ]
  },
  "sa_score": {
    "function": "step",
    "params": {
      "cutoff": 4.0,
      "scale": 4.0
    }
  }
}
```

### Key Decisions

1. **Vina 对接参数必须与比赛一致**：exhaustiveness=8, box_size=30, center=[18.28, 2.31, 21.44]
2. **Suzuki 修复已应用**：route_validity_score 预期从 0.5 → ~1.0
3. **校准提交用当前最优分子集**：不浪费提交次数
4. **归一化函数形式优先选最简单的**：避免过拟合

### File Changes

| File | Action | Description |
|------|--------|-------------|
| `src/scorer.py` | New | 本地评分引擎 |
| `src/calibrator.py` | New | 校准器 |
| `data/calibration.json` | New | 持久化校准参数 |
| `scripts/score_local.py` | New | CLI 入口 |
