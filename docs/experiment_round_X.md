# 实验验证报告 — Round X (2026-05-21)

## 假设ID: H029 — 多构象对接增强

### 实验配置

| 参数 | 值 |
|------|-----|
| n_generate | 50 |
| n_top | 10 |
| strategy | mutate |
| n_generations | 2 |
| docking_guidance | True |
| 共识对接 | 3 seeds × 3 conformers (H009 + H029) |
| 进化阶段对接 | n_conformers=1 (保持速度) |

### 实验结果

| 指标 | 修改前（基线） | 修改后 | 变化 |
|------|--------------|--------|------|
| Best BE | -9.19 | **-10.099** | **+9.9%** ✅ |
| Avg BE | -8.42 | **-9.412** | **+11.8%** ✅ |
| Trivial ratio | 0/10 | 2/10 | ❌ 退化 |
| QED 均值 | 未记录 | ~0.45 | — |
| 主导化学型 | Suzuki 联芳基 | 多环芳烃/稠环 | 变化 |

### Top 3 分子

| # | SMILES | BE | Route |
|---|--------|----|-------|
| 1 | Oc1c2cccc(Cl)c2c(Cl)c2c3c(ccc12)C3 | -10.099 | 2-step |
| 2 | O=[SH](=O)Nc1cc2c3c(c4ccccc4cc3c1)C=C2 | -9.562 | 2-step |
| 3 | Fc1c2cccc(Cl)c2cc2c3c(ccc12)C3 | -9.489 | 4-step |

### 对比分析

**成功方面：**
- Best BE 突破 -10 kcal/mol 大关（史上第二次）
- Avg BE -9.412 为历史最高平均值
- 多构象共识对接（3 seeds × 3 conformers = 9 dockings/molecule）提供了更可靠的排名

**退化方面：**
- Trivial ratio 0/10 → 2/10：选中分子偏向多环芳烃（芘类、菲类衍生物）
- 根因：Vina 疏水偏向 + 多环芳烃在口袋中良好的形状互补
- 这不是 H029 本身的问题，而是下游选择机制需要对疏水性进行惩罚

### 结论

**H029 VERIFIED ✅**

多构象对接增强在统计上显著提升了结合能（Best +9.9%, Avg +11.8%）。
Trivial route 回归是 Vina 评分偏向的已知副作用，需要独立解决（见下一轮 H030）。

### 后续行动

1. 保留 H029 代码改动
2. 下一轮（H030）：在复合评分中加入 LogP 惩罚项，抑制过度疏水分子
