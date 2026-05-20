# 实验报告 — Round 18 (H025)

## 假设: H025 — 扩展取代基库

### 实验配置
- n_generate=50, n_top=10, strategy=mutate, n_generations=2
- docking_guidance=True
- 修改: `_add_substituent()` SUBSTITUENTS 5→10

### 结果对比

| 指标 | 基线 | H025 | Δ |
|------|------|------|-----|
| Best BE | -9.153 | -8.972 | -0.181 ❌ |
| Avg BE | -8.453 | -8.325 | -0.128 ❌ |
| Trivial | 0/10 | 0/10 | 0 |

### 分析
1. 极性药效团的添加未能提升 Vina 结合能
2. Vina 评分函数偏向疏水/平面芳环体系，极性 hinge binder（BE ~ -8.05）不如萘联苯（BE -8.97）
3. 新增极性基团可能通过增加分子量/LogP导致更多分子被过滤

### 结论: ❌ REJECTED
- 不要在 Vina 对接场景下追求极性氢键优化
- 应聚焦疏水芳香体系的扩展和多样化
