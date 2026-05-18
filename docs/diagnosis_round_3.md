# 瓶颈诊断报告 — Round 3

## H012+H013 后状态
- ✅ BE: -9.755 (接近基线 -9.941, 差距仅 1.9%)
- ✅ Trivial: 0/10 
- ✅ 路线质量: 全部 ≥ 0.70
- ⚠️ Crossover 成功率仅 30%

## 剩余瓶颈

### 瓶颈 1: 种子选择仅看 BE，忽略 SA
MOOSE-Chem 强调"Multi-objective optimization: docking score AND synthetic accessibility co-optimized."
当前在种子选择阶段只用 BE，SA 仅做硬过滤。

### 瓶颈 2: Crossover 成功率可提升
30% 成功率意味着 70% 的 crossover 尝试被浪费。可通过重试机制提升。

## 假设 H014: SA 软引导 + Crossover 重试

```
假设ID: H014
瓶颈: BE种子选择无SA意识，crossover成功率低
文献支撑:
  - MOOSE-Chem (Yang 2025): Multi-objective co-optimization of BE and SA
  - MolLEO (Wang 2024b): 进化算法中重试机制提升有效变异率

───────────────── 推理链 ─────────────────
步骤 | 内容                                    | 置信度
S1   | SA软引导优先低SA分子为种子               | 高
S2   | 低SA种子产生更合成友好的后代              | 中
S3   | Crossover重试提升有效后代数               | 高

综合置信度: 中

───────────────── 验证标准 ─────────────────
Q3: 最低可接受结果？
答: BE下降不超过2%，trivial保持0/10

───────────────── 改进方案 ─────────────────
改动文件: src/generator.py
改动内容:
  1. generate_with_docking_guidance: BE种子选择时加入 SA×0.1 惩罚
  2. Crossover 重试：失败时最多重试3次
验证指标: BE + trivial + 路线质量
```
