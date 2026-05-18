# 代码演进日志 — Round 1

## 假设ID: H011

**假设**: 进化选择仅按结合能排序导致过早收敛，引入多样性保持选择可避免此问题。

## 修改文件

### `src/generator.py`

**新增函数:**
- `_diverse_selection(candidates, n_select, diversity_weight)` — 贪心 MMD 多样性选择
- `_avg_pairwise_similarity(molecules)` — 计算平均成对 Tanimoto 相似度

**修改函数:**
- `generate_with_docking_guidance()` — 种子选择阶段加入多样性保持

## 修改内容摘要

1. **种子选择逻辑改变** (原 line 634-635, 新 line ~750-790):
   - 原: `current_seeds = docked_batch[:top_k]` (纯BE选择)
   - 新: 50% 种子来自纯BE，50% 来自多样性选择
   - 多样性选择使用贪心 MMD 算法 (Maximum Minimal Distance)
   - 综合评分: `(1 - 0.5) × BE_norm + 0.5 × diversity_norm`

2. **最终分子选择增强** (新 line ~800-830):
   - 从 top 2×n_molecules 候选池中选择
   - 60% 纯BE + 40% 多样性
   - 确保最终输出既高结合能又有多样性

3. **多样性指标记录** (新 line ~795):
   - 每代打印 avg_pairwise_sim 指标
   - 日志标签: `[H002+H011]`

## 文献依据

- **MOOSE-Chem** (Yang et al., 2025): "Diverse initial population is essential for evolutionary search to avoid premature convergence"
- **MolLEO** (Wang et al., 2024b): LLM-based multi-objective evolutionary optimization
- **ChemCrow** (Bran et al., 2024): 化学空间覆盖度决定 Agent 探索边界

## 代码自检

- `python3 -m py_compile src/generator.py` ✅
- All imports verified ✅
