# 瓶颈诊断报告 — Round 10

> 基于 9 轮迭代的完整代码审查与实验数据分析，识别当前剩余瓶颈

---

## 一、当前性能基线（Round 9 / H017）

| 指标 | 当前值 | 历史最佳 | 天花板估计 |
|------|--------|---------|-----------|
| 最佳结合能 | -9.591 kcal/mol | -9.941 (H011) | ~-10.0 (Vina 精度极限) |
| Top-10 平均结合能 | -8.853 kcal/mol | -8.961 (Round 4) | ~-9.0 |
| Trivial route 比例 | 0.0 (0/10) ✅ | 0.0 (H015) | 已解决 |
| QED 均值 | ~0.55 | ~0.60 | ≥ 0.5 |
| 分子多样性 | 未系统量化 | — | **当前主要瓶颈** |

---

## 二、瓶颈诊断

### 瓶颈1: 分子生成策略多样性不足（高影响）

**现状**: generator.py 有三种生成策略：
- `mutate`: 从 SCAFFOLDS 库随机变异（主力，已验证）
- `combine`: **直接字符串拼接 SMILES**（如 `c1ccc(cc1)C(=O)Nc1ccccc1`），仅在巧合时才产生有效分子
- `random`: 等同于高强度 mutate

**根本原因**:
当前 `combine` 策略的实现 (generator.py:664-685)：
```python
combined = "".join(parts)  # ← 字符串拼接，非化学操作
mol = Chem.MolFromSmiles(combined)
```
这不是真正的分子重组——它只能拼接预定义骨架和连接子的 SMILES 字符串。对于复杂片段（如含环闭合、多连接点），必然失败。更重要的是，它**只使用硬编码的 SCAFFOLDS 列表**，无法利用对接成功分子的特征片段。

**为什么重要**:
MOOSE-Chem 明确指出："Diverse initial population is essential for evolutionary search to avoid premature convergence." 当前 docked 分子池中的高分分子，其 BRICS 片段代表已验证的高质量化学型。用这些片段重组可以比随机变异更高效地探索相关化学空间。

**文献支撑**:
- MOOSE-Chem (Yang et al., 2025): 进化算法中重组(recombination)算子是产生多样性的核心
- MolLEO (Wang et al., 2024b): LLM 驱动的重组操作，多种进化算子组合
- Deep Lead Optimization (JACS 2024): BRICS 分解→重组是 Fragment Replacement 的基石
- ChemCrow (Bran et al., 2024): 工具/策略多样性决定 Agent 能力边界

### 瓶颈2: Crossover 算子利用率低（中影响）

**现状**: 在 docking guidance 模式中，crossover 仅以 15% 概率使用 (generator.py:846)。MOOSE-Chem 和 MolLEO 都强调重组算子的关键作用。

**根本原因**: 初始设计时变异作为主力算子，crossover 是后来加入的（H013），概率保守。

**文献支撑**: MOOSE-Chem: "crossover between parent molecules, fragment swapping, and scaffold hopping"

### 瓶颈3: 缺失 Fragment Replacement 算子（中影响）

**现状**: Deep Lead Optimization 定义四个核心子任务中，Fragment Replacement 完全未实现。当前只有 Scaffold Hopping (H010)。

**文献支撑**: JACS 2024: Fragment replacement 是先导化合物优化四大核心操作之一

---

## 三、假设提出

### H018: 修复并增强 BRICS 片段重组 "combine" 策略

```
假设ID: H018
瓶颈: generator.py "combine" 策略是字符串拼接，不产生真正的化学多样性。
      更重要的是，它只从硬编码 SCAFFOLDS 库取片段，
      完全未利用对接成功分子（已验证的高质量化学型）中的 BRICS 片段。
文献支撑:
  - MOOSE-Chem (Yang et al., 2025): "Diverse initial population is essential"
  - MolLEO (Wang et al., 2024b): 重组算子是进化的核心驱动力
  - JACS 2024: BRICS 分解→重组是 Fragment Replacement 的标准方法

───────────────── 推理链 ─────────────────
步骤 | 内容                                            | 置信度 | 推理方式 | 依据来源
S1   | 当前 combine 策略只生成有限的有效分子               | 高     | 代码审查 | generator.py:664-685
S2   | BRICS Build 可以从片段库生成有效新分子              | 高     | 实验验证 | BRICSBuild 返回4个有效分子
S3   | 从已对接分子提取 BRICS 片段 → 高质量化学型种子       | 高     | 演绎     | 对接筛选 = 化学型已验证
S4   | BRICS 重组产生的新分子可能具有 >= 原有结合能         | 中     | 演绎     | 片段来自高分分子但组合可能不兼容
S5   | 增加化学多样性 → 提高找到更优分子的概率             | 中     | 演绎     | MOOSE-Chem 多样性原则

综合置信度 = 中

───────────────── 验证标准 ─────────────────
Q1 如果结合能未提升但多样性增加，是否仍保留？
答：是
理由：更高多样性是长期进化搜索的基础，即使单轮不提升结合能

Q2 如果指标下降，最可能的原因是什么？
答：BRICS Build 生成大量无效/低质量分子（BRICS Build 不加约束会组合爆炸）
理由：BRICS Build 不检查类药性质，需后续过滤

Q3 本假设的最低可接受结果是什么？
答：top-10 平均结合能下降不超过 3%，且生成池有效分子率 ≥ 20%

───────────────── 改进方案 ─────────────────
改动文件: src/generator.py
改动内容:
  1. 新增函数 _brics_recombine(pool_smiles_list, n_molecules):
     - 从已有分子池做 BRICS 分解，收集片段
     - 随机采样片段组合，用 BRICS Build 生成新分子
     - 过滤：有效性 + 类药性质 + 去重
  2. 修复 combine 策略：不再字符串拼接，改为使用 BRICS 片段库
  3. 在 generate_with_docking_guidance 中：每次从当前 docked 池
     提取片段，为 combine 策略提供"活的"片段库
  4. 将 crossover 使用率从 15% 提升到 25%

验证指标: top10 平均结合能、QED 分布、有效分子率、分子多样性 (avg pairwise similarity)
```

---

## 四、其他待验证假设

### H019: 提升 Crossover 算子概率（已合入 H018）
- 15% → 25%，简单参数调整

### H020: Fragment Replacement 算子
- 影响: 中高（补全 JACS 四大核心操作）
- 难度: 中（需新写 ~100 行）
- 优先级: Round 11

---

## 五、本轮行动计划

1. **阶段三（代码演进）**: 实施 H018 — BRICS 片段重组策略
2. **阶段四（实验验证）**: 运行 run_pipeline (n_generate=50, n_top=10) 验证
3. **复盘**: 分析新策略的分子多样性、有效率和结合能变化
