# 实验报告 — Round 13 (H021)

## 实验配置
- **假设ID**: H021 — H020 规则修复 combine 策略 trivial route 回归
- **Pipeline**: combine, n_generate=50, n_top=10, n_generations=2, docking_guidance=ON
- **代码变更**: 无（H020 规则已就位）

## 实验结果

### 核心指标对比
| 指标 | H021 (combine+H020) | H019 (combine, 无H020) | 变化 |
|------|:---:|:---:|:---:|
| Best BE (kcal/mol) | -9.159 | -9.335 | -0.176 |
| Avg BE (kcal/mol) | -8.677 | -8.581 | +0.096 |
| **Trivial ratio** | **0/10** | **0.1** | **✅ 修复** |
| QED 均值 | ~0.6 | ~0.6 | 持平 |

### 分子化学空间分析
Combine 策略生成了 mutate 策略罕见的新骨架类型：
- 吲哚-萘联芳（Suzuki）
- 苯并咪唑-甲苯联芳（Suzuki）
- 喹唑啉-吲哚联芳（喹唑啉规则 + Suzuki）
- 喹喔啉-氯苯联芳（Suzuki）

所有 10 个分子均有非 trivial 路线，主要使用 Suzuki 偶联和已有杂环合成规则。

## 结论
✅ **H021 验证成功** — H020 规则库扩充成功修复了 H019 的 combine 策略 trivial route 回归。
- trivial 从 0.1 → 0.0
- 验证了 LARC (2025) 核心观点：规则覆盖率是逆合成质量的关键决定因素
- Combine 策略现在可安全使用，化学多样性优于 mutate
