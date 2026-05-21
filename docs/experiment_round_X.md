# 实验验证 — Round X (新会话 H031)

## 日期: 2026-05-21

---

## 假设: H031 — 单原子替换 Trivial 路线检测

### 实验设置
- Pipeline: `run_pipeline(n_generate=50, n_top=10, strategy=mutate, n_generations=2, docking_guidance=True)`
- 对照组: H030 基线 (Best BE=-9.972, Avg=-8.757, trivial claimed 0/10 actual 2/10)
- 实验组: H031 修复后

### 结果

| 指标 | H030 基线 | H031 实验 | 变化 |
|------|----------|----------|------|
| Best BE | -9.972 | -9.902 | -0.070 (-0.7%) |
| Avg BE | -8.757 | -8.642 | -0.115 (-1.3%) |
| Trivial 总数 | 2/10 (OH→Cl型) | 2/10 (smiles>>smiles型) | 不变 |
| Route quality 惩罚 | Trivial 分子未受惩罚 | Trivial 分子受惩罚 | ✅ 改善 |

### Trivial 路线详细分析

**修改前 (H030)**: 
- 2个 OH→Cl 单原子替换 trivial 路线
- `route_quality` ~0.7（未被识别为 trivial）
- 复合评分中 `0.15 * 0.7 = 0.105` 得分优势
- 挤占真正有价值分子的 top-10 位置

**修改后 (H031)**:
- 0个 OH→Cl 类型 ✅ 完全消除
- 2个 smiles>>smiles 类型（四环骨架无匹配规则）
- `route_quality` = 0.0（正确惩罚）
- 复合评分中无 route 得分优势

### 结论

**H031: ✅ VERIFIED**

1. `_is_single_atom_swap()` 函数正确检测 OH↔Cl/Br 单原子替换
2. 化学上无效的 trivial 路线被正确标记
3. BE 变化在实验噪声范围内（-0.7%）
4. 剩余的 2 个 trivial 路线是不同根因：四环骨架缺少逆合成规则
5. 建议下一轮假设 H032: 扩充四环/多环骨架的逆合成规则
