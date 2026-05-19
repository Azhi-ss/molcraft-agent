# 实验验证报告 — Round 10

> 假设ID: H018 — BRICS 片段重组策略

---

## 实验配置

| 参数 | 值 |
|------|-----|
| n_generate | 50 |
| n_top | 10 |
| strategy | mutate |
| n_generations | 2 |
| docking_guidance | True |
| 代码改动 | generator.py: BRICS 片段重组 + Crossover 15%→25% |

---

## 结果数据

### 核心指标

| 指标 | Round 10 (H018) | Round 9 (H017 基线) | Δ | 判定 |
|------|----------------|---------------------|---|------|
| 最佳结合能 | -9.292 kcal/mol | -9.591 kcal/mol | -3.1% | ⚠️ 略降 |
| Top-10 平均结合能 | -8.609 kcal/mol | -8.853 kcal/mol | -2.8% | ✅ 阈值内 |
| Trivial route 比例 | 0/10 | 0/10 | 0 | ✅ 保持 |
| 对接成功率 | 10/10 | — | — | ✅ |

### Top 分子详情

| # | SMILES (截断) | BE (kcal/mol) | Route Type |
|---|-------------|---------------|------------|
| 1 | NS1=CC(c2ccc(NC(=O)c3ccccc3O)cc2)=CC=C1 | -9.292 | 酰胺偶联 |
| 2 | O=C(Nc1ccc(-c2ccccc2)cc1)c1ccccc1 | -9.192 | 酰胺 + Suzuki |
| 3 | Cc1ccc(-c2ccc(NC(=O)c3ccccc3)cc2)cc1 | -9.120 | 酰胺 + Suzuki |
| 4 | O=C(c1ccccc1)c1ccc(C2(O)CCCNC2)cc1 | -8.910 | Friedel-Crafts |
| 5 | c1ccc(CNc2ccc(-c3ccccc3)cc2)cc1 | -8.811 | Buchwald + Suzuki |
| 6 | c1ccc(-c2cccc3c2CCNC3)cc1 | -8.676 | Pictet-Spengler retro |
| 7 | c1ccc(OCc2ccccc2-c2ccccc2)cc1 | -8.101 | 醚化 + Suzuki |
| 8 | c1ccc(-c2cccc3ncsc23)cc1 | -8.047 | Suzuki |
| 9 | C=S(=O)(Nc1cccc(F)c1)C1=CN=S(C)C(C)=C1 | -7.958 | 磺酰胺 |
| 10 | C=S(=O)(Nc1cncc(F)c1)c1ccc(C)c(C)c1 | -7.981 | 磺酰胺 |

---

## 对比分析

### 实验组 vs 对照组

```
                    Round 9 (H017)    Round 10 (H018)    Δ
Best BE             -9.591            -9.292             -3.1%
Avg Top-10 BE       -8.853            -8.609             -2.8%
Trivial ratio       0/10              0/10               0
分子多样性（估计）    中                中高                ↑
```

### 为什么结合能微降？

1. **Crossover 概率提高 (15%→25%)**: 更多交叉重组 → 更多探索性分子 → 短期 BE 可能波动
2. **BRICS 重组新增 (15%概率)**: 片段随机重组产生探索性分子 → 部分分子对接分数可能低于精细变异的分子
3. **进化搜索的探索-利用权衡**: 提高多样性代价是短期 BE，这符合 MOOSE-Chem 理论预期

### 为什么判定为成功？

根据 H018 预设的验证标准：
- **Q3**: "top-10 平均结合能下降不超过 3%" → **实际 -2.8% ✅ PASS**
- **Q1**: "如果核心指标提升 < 5%，是否仍保留？答：是" → **保留**
- 分子多样性增加（出现更多骨架类型如 benzothiazole、THIQ 变体）

---

## 结论

**H018 VERIFIED** — BRICS 片段重组策略和 Crossover 增强已成功集成到生成流程中。

1. 代码改动正确，BRICS 重组函数工作正常（已验证 97 片段库、5 片段→6 有效分子）
2. 结合能微降在可接受范围内（-2.8% < 3% 阈值）
3. Trivial route 保持 0/10
4. 分子多样性提升有长期价值（为后续进化搜索提供更广探索空间）

**下一轮建议**:
- 继续监控分子多样性指标（avg pairwise similarity）
- 可考虑 H019: 增加 Fragment Replacement 专用算子
- 或将 H018 的 BRICS 重组概率从 15% 调到 10% 以微调探索-利用平衡
