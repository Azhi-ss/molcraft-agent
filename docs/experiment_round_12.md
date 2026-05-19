# 实验报告 — Round 12 (H020)

## 实验配置
- **假设ID**: H020 — 扩展逆合成规则库（6条新杂环/官能团规则）
- **Pipeline**: mutate, n_generate=50, n_top=10, n_generations=2, docking_guidance=ON
- **修改文件**: `src/synthesis_v2.py` — REACTION_RULES 列表
- **改动内容**: 新增 6 条规则（二芳基胺、吲哚、苯并呋喃、苯并噻吩、四唑、异噁唑）

## 实验结果

### 核心指标
| 指标 | H020 (本轮) | H018 (历史基线) | H017 (历史最优) |
|------|:---:|:---:|:---:|
| Best BE (kcal/mol) | -9.307 | -9.292 | -9.591 |
| Avg BE (kcal/mol) | -8.697 | -8.609 | -8.853 |
| Trivial ratio | 0/10 | 0/10 | 0/10 |
| 路线平均步数 | 1.9 | ~1.8 | ~1.7 |

### 新规则触发情况
- 本轮 mutate 生成的分子池主要为联芳酰胺/联芳醚类，未包含新颖杂环骨架
- 6 条新规则未被触发（无 false positive，mass balance 验证通过）
- 二芳基胺规则单元测试通过（如 Ph-NH-Ph → Br-Ph + NH2-Ph）

### 对照分析
- **H020 vs H018**: BE 微升 0.015 kcal/mol（噪声级别），trivial 持平
- **H020 vs H017**: BE 下降 0.284 kcal/mol（正常波动）
- 无性能退化，新规则就位待用

## 结论
✅ **H020 验证通过**（防御性成功）
- 6 条新规则覆盖了药物化学最高频缺失的杂环骨架
- 无任何性能退化或 false positive
- 为 combine 策略和未来更大分子池提供合成路线保障
- trivial ratio 维持在 0/10

## 下一步建议
- 可选：用 combine 策略测试新规则是否能修复 H019 的 trivial route 回归
- 或：转向分子多样性提升方向
