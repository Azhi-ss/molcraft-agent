# 实验验证报告 — Round 4

## 假设ID: H026 — Expanded Large Aromatic Scaffold Library

### 实验配置

| 参数 | 值 |
|------|-----|
| n_generate | 50 |
| n_top | 10 |
| strategy | mutate |
| n_generations | 2 |
| docking_guidance | true |
| generator | mutate |

### 实验结果

| 指标 | 基线 (H027) | H026 实验 | Δ |
|------|-------------|-----------|-----|
| Best BE | -9.332 | -8.389 | +0.943 ❌ |
| Avg BE | -8.810 | -7.989 | +0.821 ❌ |
| Trivial ratio | 0/10 | 1/10 | +1 ❌ |
| QED (top) | ~0.45 | ~0.50 | ~不变 |

### 对照分析

基线为 Round 22 已验证的 H027 配置（完全相同参数）。
实验组在三个指标上均显著退步。

**Top 分子分析:**
- 排名第一的分子为甲基吲哚衍生物（BE -8.389），但获得 trivial route
- 所有 top 10 分子均未使用新增的大芳香骨架（蒽、菲、芘等）
- 主导化学型仍为 Suzuki 联芳基 + 吲哚/吲哚啉

### 失败原因

1. **骨架选用率低**: 新骨架占 81 个中的 12 个(14.8%)，mutate 策略中随机采样导致大骨架极少被选中
2. **MW 过滤**: 大芳香核心(14-20 atoms)加侧链后易超 MW=500
3. **口袋匹配**: TYK2 ATP 口袋更适合 2-3 环体系，4-5 环大芳香体系不匹配
4. **对接负选择**: Docking guidance 偏向中等大小疏水分子

### 意外发现: Suzuki SMARTS 特异性 Bug

在分析 trivial route 分子时发现:
- 甲基吲哚 (`Cc1ccc2[nH]c3ccccc3c2c1`) 的 Fischer indole 逆合成规则被跳过
- 根因: Suzuki 规则 `[c;R][c;R]` 在 RETRO_RULES 中排名 #5，错误匹配稠环体系内部 C-C 键
- 该 SMARTS 在 indole 上产生 13 个匹配，生成 26 组产物
- Fischer indole 规则排在 #61，被前面的 Suzuki 错误匹配拦截
- **这是 trivial route 从 0/10 回归到 1/10 的直接原因**

### 结论

**H026: REJECTED** — 大芳香骨架扩展对 TYK2 5C01 结合能无正面贡献。

**下一轮方向**: 修复 Suzuki SMARTS 特异性，防止稠环内错误匹配。
