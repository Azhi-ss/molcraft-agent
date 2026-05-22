# 诊断报告 — Round 3 (H036)

## 瓶颈诊断

多构象对接（H035）已显著提升 Avg BE (+7.2%), 但对接引导的探索-利用平衡可进一步优化。

当前参数: batch_size=10, top_k=5 (50% 保留率), n_generations_inner=3
问题: top_k=5 过于激进，可能过早收敛到局部最优

## 假设 H036: 增强对接引导探索能力

增大 batch_size (10→15) + top_k (5→8) + inner generations (3→4)

预期: 更多样化的种子池 → 更优的进化收敛 → BE 再提升 0.1-0.3 kcal/mol
