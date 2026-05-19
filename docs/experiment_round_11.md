# 实验验证报告 — Round 11

> 假设ID: H019 — strategy="combine" 主模式验证

---

## 实验配置

| 参数 | 值 |
|------|-----|
| n_generate | 50 |
| n_top | 10 |
| strategy | **combine** (vs mutate in R10) |
| n_generations | 2 |
| docking_guidance | True |

---

## 结果数据

### 核心指标对比

| 指标 | Round 11 (combine) | Round 10 (mutate) | Δ | 判定 |
|------|-------------------|-------------------|-----|------|
| 最佳结合能 | -9.335 kcal/mol | -9.292 kcal/mol | +0.4% | ✅ 持平 |
| Top-10 平均结合能 | -8.581 kcal/mol | -8.609 kcal/mol | -0.3% | ✅ 阈值内 |
| Trivial route 比例 | 1/10 | 0/10 | **+1** | ❌ 退化 |
| 对接成功率 | 10/10 | 10/10 | 0 | ✅ |

### Top 分子详情

| # | SMILES (截断) | BE (kcal/mol) | Route Type | Trivial? |
|---|-------------|---------------|------------|-----------|
| 1 | O=C(Nc1cccc(-c2ccccc2)c1)c1ccccc1 | -9.335 | 酰胺 + Suzuki | No |
| 2 | Cc1ccc2c(c1N)Cc1ccccc1C2 | -9.192 | — | **YES** ⚠️ |
| 3 | O=C(Nc1ccccc1)c1ccccc1-c1ccccc1 | -9.036 | 酰胺 + Suzuki | No |
| 4 | Oc1cccc(-c2ccc3ccccc3c2)c1 | -8.898 | Suzuki | No |
| 5 | Cc1ccc(S(=O)(=O)Nc2ccccc2C2(F)CCCCC2)cc1 | -8.615 | 磺酰胺 | No |
| 6 | Cc1ccc(-c2ccc(N)cc2Cl)cc1O | -8.357 | Suzuki | No |
| 7 | O=C(c1ccccc1)c1ccccc1-c1ccccc1 | -8.212 | FC酰化 + Suzuki | No |
| 8 | O=Cc1ccc(-c2cc(O)cc(F)c2)cc1 | -8.108 | Suzuki | No |
| 9 | O=C(NC1(C2(F)COCCC2O)CCCCC1)c1ccccc1 | -8.086 | 酰胺 | No |
| 10 | c1cc2c(c(-c3cccc4c3OCCO4)c1)OCCN2 | -7.966 | Suzuki | No |

---

## 失败分析

### 为什么出现了 trivial route？

分子 #2 `Cc1ccc2c(c1N)Cc1ccccc1C2` 是一个**二苯并氮杂环庚烷**（dibenzazepine）骨架：
- 两个苯环通过一个七元含氮环桥连
- BE = -9.192 kcal/mol — **高结合能**，说明该骨架与靶点结合良好
- 该骨架由 BRICS 片段重组（combine 策略）创造 — 是一个全新的化学型

**根本原因**: 现有 35+ 条逆合成规则不覆盖二苯并氮杂环庚烷骨架的断键。该骨架的合成通常需要：
- 逆 Bischler-Napieralski 环化（已有 benzazepine 规则 H015，但苯并氮杂环庚烷的具体取代模式未覆盖）
- 或 Ullmann 偶联 + 分子内还原胺化

### BRICS 重组确实创造了新骨架

这正是 combine 策略的价值：它产出了 mutate 策略从未见过的化学型。缺点是规则库未能覆盖。

---

## 结论

**H019 REJECTED** — trivial route 标准未通过（1/10 > 0/10 阈值）。

但实验仍有科学价值：
1. **BE 性能持平**: combine 策略与 mutate 策略在结合能上表现相当（最佳 -9.335 vs -9.292, 平均 -8.581 vs -8.609）
2. **新骨架产出**: combine 策略成功创造了逆合成规则库未覆盖的新化学型（二苯并氮杂环庚烷）
3. **规则库缺口定位**: 精确识别了需要新增的规则方向（二苯并氮杂环庚烷断键）
