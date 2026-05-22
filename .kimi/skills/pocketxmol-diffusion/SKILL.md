---
name: pocketxmol-diffusion
description: PocketXMol diffusion model for pocket-aware de novo molecule generation. GPU-accelerated structure-based drug design — generates 3D molecules inside protein binding pockets. Use when you need high-quality initial seeds that already fit the pocket shape.
---

# PocketXMol Diffusion Model Usage Guide

## What It Is

PocketXMol 是一个原子级生成基础模型，通过去噪过程在蛋白口袋内生成 3D 分子结构。它直接感知口袋的几何和化学特征，生成的分子天然适配结合位点，避免了 RDKit 随机变异"大海捞针"的问题。

**已知特性（TYK2 实测）：**
- 推理时间：~25-30 秒/批次（RTX 4090）
- 生成分子类型：口袋感知，对于激酶靶点通常产出嘌呤/嘧啶核苷酸类似物
- QED 范围：0.15-0.55（含磷酸基团分子 QED 偏低是正常的）
- CFD 置信度：1.0-1.5（越高越好）

## When to Use

### ✅ 推荐使用场景

| 场景 | 模式 | 理由 |
|------|------|------|
| 新靶点第一轮探索 | `hybrid` | 扩散提供口袋感知种子 + RDKit 补充多样性 |
| 结合能卡住无提升 | `hybrid` | 扩散模型探索 RDKit 变异不了的化学空间 |
| 需要激酶铰链结合分子 | `diffusion` | PocketXMol 天然生成 hinge-binding 杂环 |
| 快速验证口袋坐标是否正确 | `diffusion` | 扩散模型生成结果可验证口袋位置合理性 |

### ❌ 不推荐使用场景

- 已有高分分子需要微调 → 用 `generate_molecules(strategy="mutate", scaffold=...)`
- GPU 服务不可用 → 自动降级无意义，等待服务恢复或用 RDKit
- 生成数量 > 100 → 扩散模型单次上限约 100 个分子

## Available Tools

### 1. `diffusion_generate` — 独立调用（两阶段优化，默认开启）

```
diffusion_generate(pdb_path="data/target.pdb", n_molecules=20, two_stage=True)
```

**两阶段流程（默认）：**
1. sbdd 生成 N 个口袋感知分子
2. 按 QED 排序，取 top-3 种子
3. 每个种子用 opt_mol 生成 5 个优化变体
4. 合并去重，返回 N + 15 个分子

**提效：** opt_mol 优化可提升 QED ~30%（实测 0.35 → 0.46），同时改善 logP 和 SA。

**参数：**
- `n_molecules`: 建议 10-50，过大无益（口袋化学空间有限）
- `two_stage`: 是否启用两阶段优化，默认 True。单阶段用 `two_stage=False`
- `n_optimize`: 两阶段时优化的种子数，默认 3
- `pocket_center`: 留空自动从 `src/config.py` 读取

### 2. `run_pipeline(generator="hybrid")` — Pipeline 集成

```
run_pipeline(n_generate=50, generator="hybrid", strategy="mutate")
```

三种模式：
- `"mutate"` — 纯 RDKit 变异（默认，不依赖 GPU）
- `"diffusion"` — 纯 PocketXMol 生成（需 GPU 服务在线）
- `"hybrid"` — 扩散种子 + RDKit 多样性补充（推荐）

## Hybrid Mode Workflow

```
扩散模型生成 n 个口袋感知分子
        +
RDKit 变异生成 m 个多样性分子
        ↓
合并候选池 → 对接 → 进化 → 逆合成 → Top N
```

扩散提供"高起点"（口袋适配），RDKit 提供"广度"（化学空间覆盖）。如果扩散模型不可用，hybrid 模式会报错而不是静默降级——这是设计意图，确保你不会在不知道的情况下跑了一个退化版本。

## GPU Server Status Check

在调用扩散模型前，可以先检查服务健康状态：
```
# 通过工具间接检查（推荐）
diffusion_generate(n_molecules=2)  # 小批量快速验证

# 或通过 Shell 直接检查
Shell: curl -s http://localhost:8001/health
```

健康响应示例：
```json
{"status":"ok", "model_loaded":true, "gpu_available":true, "active_jobs":0}
```

## Troubleshooting

| 症状 | 可能原因 | 处理 |
|------|---------|------|
| `ToolError: GPU 服务器不可达` | SSH 隧道断开 | 检查 `curl localhost:8001/health`，重建隧道 |
| `ToolError: 扩散模型推理失败` | PDB 格式问题或口袋坐标错误 | 先跑 `identify_target()` 验证口袋坐标 |
| 生成分子全是核苷酸类似物 | 正常行为——TYK2 是激酶，ATP 口袋 | 用 `hybrid` 模式补充 RDKit 多样性 |
| 推理超时（>30min） | GPU 服务器过载或 PDB 太大 | 减小 n_molecules，检查服务器 `nvidia-smi` |
| 生成分子 QED 很低（<0.2） | 含磷酸基团分子 QED 天然低 | 这些分子对接可能很好，不要仅凭 QED 过滤 |

## Interaction with Other Skills

- **molcraft-vina-strategies**: 扩散模型生成的分子对接时，口袋坐标必须一致。如果 `identify_target()` 建议调整 `DOCKING_CENTER`，扩散模型会自动使用新坐标。
- **molcraft-synthesis-rules**: 扩散模型生成的分子可能含磷酸酯/核苷结构，这些结构在 `synthesis_v2.py` 中可能缺少 SMARTS 规则，需要扩充。
- **molcraft-hypothesis-template**: 扩散模型引入是一个有效假设，验证指标应该包含"扩散组 vs RDKit 组的 top10 平均结合能差值"。
