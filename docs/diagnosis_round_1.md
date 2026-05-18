# 瓶颈诊断与假设报告 — Round 1

## 1. 现有基线性能

| 指标 | Round 1 (无 DG) | Round 2 (有 DG) | Round 3 (有 DG) |
|------|:---:|:---:|:---:|
| 最佳结合能 | -8.117 | -8.430 | -8.335 |
| Top-10 平均结合能 | -7.888 | -7.901 | -8.131 |
| Trivial route 比例 | 0% | 0% | 0% |
| Docking guidance | OFF | ON | ON |

**基线结论**: 结合能稳定在 -8.1 ~ -8.4 kcal/mol，合成路线 100% 有效（0% trivial）。核心可改进维度是**结合能**。

## 2. 瓶颈诊断

### 瓶颈 1（主要）: 进化选择仅按结合能排序 → 过早收敛

**根因分析**:
- `generate_with_docking_guidance()` (generator.py:627-635) 每代选择种子时：
  ```python
  docked_batch.sort(key=lambda x: x.get("binding_energy", 999))
  current_seeds = docked_batch[:top_k]
  ```
- 纯粹按结合能选 top_k 种子，不考虑分子间相似度
- 后果：经过 1-2 代后，所有种子分子高度相似（Tanimoto > 0.7），进化停滞于局部最优

**文献支撑**:
- **MOOSE-Chem** (Yang et al., 2025): "Diverse initial population is essential for evolutionary search to avoid premature convergence"
- **MolLEO** (Wang et al., 2024b): 多目标进化优化平衡 exploration vs exploitation
- **ChemCrow** (Bran et al., 2024): 化学空间覆盖度决定探索边界

### 瓶颈 2（次要）: 分子生成多样性不足

**根因分析**:
- 变异算子虽然是随机的，但都是从固定库中选取取代基/骨架
- 缺少系统性的侧链装饰（Side-Chain Decoration），Deep Lead Optimization 定义为核心子任务

### 瓶颈 3（低优先级）: 逆合成规则覆盖率

**根因分析**: trivial route 比例已为 0%，说明当前 35 条规则足够。提升空间有限。

## 3. 假设 H011: 多样性保持选择

```
假设ID: H011
瓶颈: 进化选择仅按结合能排序，导致分子过早收敛于局部最优
文献支撑:
  - MOOSE-Chem (Yang et al., 2025): "Diverse initial population is essential
    for evolutionary search to avoid premature convergence"
  - MolLEO (Wang et al., 2024b): 多目标进化优化中多样性是核心维度
  - ChemCrow (Bran et al., 2024): 化学空间覆盖度 = Agent 探索边界

───────────────── 推理链（拆解假设的推理过程） ─────────────────
步骤 | 内容                                          | 置信度 | 推理方式 | 依据来源
S1   | 仅按结合能选择种子导致分子高度相似               | 高     | 演绎     | 进化算法理论
S2   | 分子高度相似后变异难以跳出局部最优               | 高     | 文献     | MOOSE-Chem 2025
S3   | 引入多样性奖励（Tanimoto 距离）可保持探索性       | 高     | 文献     | MolLEO 2024b
S4   | 保持探索性可发现结合能更优的新化学空间             | 中     | 演绎     | 从 S3
S5   | 最终 top-10 平均结合能提升                       | 中     | 演绎     | 从 S4

综合置信度 = 所有步骤中最低的（木桶效应）：中
S4/S5 为"中"——探索更广空间不一定保证找到更好的分子，可能只是发现更多中等质量分子

───────────────── 验证标准（提前定好成败门槛） ─────────────────
Q1 如果核心指标提升 < 5%，是否仍保留？
答：是
理由：即使结合能不变，分子多样性提升本身有价值（为后续优化提供更广种子池）

Q2 如果指标下降，最可能的原因是什么？
答：多样性奖励权重过大，导致低结合能分子占据种子位置，降低选择压力
理由：exploration-exploitation tradeoff 是进化算法的经典挑战

Q3 本假设的最低可接受结果是什么？
答：top-10 平均结合能下降不超过 3%（即不低于 -7.887），且分子间平均 Tanimoto < 0.5（当前预计 > 0.7）

───────────────── 改进方案 ─────────────────
改动文件: src/generator.py
改动内容:
  1. 在 generate_with_docking_guidance 的种子选择阶段，引入多样性保持算法:
     - 使用 Maximum Minimal Distance (MMD) 或 Determinantal Point Process (DPP) 选择
     - 实现方法: 贪心选择 — 先选结合能最优分子，然后迭代选择与已选集合
       Tanimoto 距离最远的分子（结合能作为 tie-breaker）
  2. 同时保留 top_k 个纯结合能最优分子 + (n_seeds - top_k) 个多样性分子
  3. 最终 top-N 选择也加入多样性考量
验证指标:
  - 核心: top-10 平均结合能（对比基线 -8.131）
  - 辅助: 分子间平均 Pairwise Tanimoto 距离、QED 均值、trivial route 比例
  - 运行 pipeline: n_generate=50, n_top=10, strategy=mutate, docking_guidance=ON
```
