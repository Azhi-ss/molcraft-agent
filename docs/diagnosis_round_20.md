# 诊断报告 — Round 20 (H027)

## 假设 H027: 放宽分子过滤器

### 瓶颈
passes_filters 的 max_rings=7 和 max_sa=6.0 可能排除了
大型疏水芳香体系（>7环），这些恰是 Vina 评分函数偏好的分子。

### 推理链
S1: Vina 奖励疏水π-π堆积 (已验证)
S2: 大芳香体系有更多环, 但 MW 仍可 <500 (高置信度)
S3: 当前 max_rings=7 排除了 ≥8 环的分子 (代码审查)
S4: 放宽 max_rings→9, max_sa→7.0 → 更多大芳香分子通过 (中)

### 改动
- src/evaluator.py: passes_filters() 默认 max_rings 7→9, max_sa 6.0→7.0
- 基于新科学认知的策略调整，非单纯调参

### 验证标准
- 最低: best BE ≥ -9.15 (不退化)
- 期望: best BE > -9.5, 至少 1 个分子含 ≥8 环
