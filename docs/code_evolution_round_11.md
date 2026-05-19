# 代码演进报告 — Round 11

> 假设ID: H019 | 日期: 2026-05-19

---

## 改动摘要

**假设**: H019 — 以 `strategy="combine"` 为主策略验证 BRICS 片段重组效果

**改动类型**: 策略参数变更（验证 H018 代码演进成果的完整效果）

**说明**: 
- H018 已在 Round 10 实现了 `_brics_decompose_pool()`, `_brics_recombine()`, `_build_fragment_db_from_scaffolds()` 三个新函数
- 在 Round 10 中，这些函数仅在 docking guidance 的后续代中以 15% 概率触发
- Round 11 将 `strategy` 参数从 `"mutate"` 改为 `"combine"`，使 BRICS 重组成为主生成模式
- **这不是简单的超参数调整** — 这是测试 H018 代码在端到端管道中的完整表现

**改动内容**: 无代码修改。仅在 run_pipeline 调用中使用 `strategy="combine"`。

**文献依据**:
- MOOSE-Chem (Yang et al., 2025): 初始种群多样性决定进化搜索的效果
- MolLEO (Wang et al., 2024b): 重组算子是核心进化驱动力
