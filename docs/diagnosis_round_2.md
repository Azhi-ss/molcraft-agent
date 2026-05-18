# 瓶颈诊断报告 — Round 2

## H012 后状态
- ✅ 路线质量: 不现实吡啶路线 3个→0个
- ✅ 路线评分复合选择正常工作
- ⚠️ BE: -9.335 (基线 -9.941, -6.1%)
- ⚠️ 仅 5 种变异算子，缺少双亲重组(crossover)

## 新瓶颈: 分子探索空间受限

**现象**: 当前所有变异算子为单亲操作（add/replace/remove/insert/scaffold_hop），
无法组合两个高分分子的有益片段。这是进化算法的已知限制。

**根因**: 缺少 crossover（交叉重组）算子，限制了化学空间的有效探索。

**文献**: MOOSE-Chem (Yang 2025) — "crossover between parent molecules" 是核心变异算子；
MolLEO (Wang 2024b) — LLM 驱动的重组操作能显著提升分子多样性。

## 假设 H013: 分子 Crossover 重组算子

```
假设ID: H013
瓶颈: generator.py 只有单亲变异，缺少双亲交叉重组
文献支撑:
  - MOOSE-Chem (Yang et al., 2025): "Evolutionary operators include crossover
    between parent molecules, fragment swapping, and scaffold hopping"
  - MolLEO (Wang et al., 2024b): LLM 驱动的重组操作提升化学空间探索效率
  - Deep Lead Optimization (JACS 2024): Fragment replacement 是先导优化核心子任务

───────────────── 推理链 ─────────────────
步骤 | 内容                                    | 置信度 | 推理方式 | 依据来源
S1   | Crossover 引入新化学空间区域            | 高     | 文献     | MOOSE-Chem
S2   | 新空间区域可能含更高 BE 分子             | 中     | 演绎     | 进化算法理论
S3   | 双亲重组不会破坏现有优秀分子             | 高     | 演绎     | 仅影响新生成分子

综合置信度: 中

───────────────── 验证标准 ─────────────────
Q1: 核心指标提升 < 5%，是否仍保留？
答: 否
理由: H013 目标是 BE 提升。若 BE 无明显改善，crossover 算子价值有限

Q2: 指标下降最可能原因？
答: Crossover 产生无效分子过多，浪费对接资源
理由: 随机交换片段可能产生化学上不合理的结构

Q3: 最低可接受结果？
答: 最佳 BE ≥ -9.5 或平均 BE 提升 ≥ 0.2 kcal/mol，trivial 不增加

───────────────── 改进方案 ─────────────────
改动文件: src/generator.py
改动内容:
  1. 新增 _crossover_mol(mol1, mol2) 函数
  2. 在 _mutate_mol 中增加 ~15% 概率调用 crossover
  3. Crossover 算法: BRICS 分解两分子 → 交换片段 → 重组验证
验证指标: BE + 多样性(avg_pairwise_sim) + trivial 比例
```
