# 代码演进报告 — Round 4

## 假设ID: H026

**假设**: 扩展大芳香环骨架库（新增12个3-5环多环芳烃/杂芳烃骨架）

### 修改内容

**文件**: `src/generator.py`
**位置**: `SCAFFOLDS` 列表（第38-136行 → 扩展至第171行）
**改动类型**: 纯增量（只添加，不修改现有逻辑）

### 新增骨架明细

| # | 名称 | 环数 | 原子数 | SMILES |
|---|------|------|--------|--------|
| 1 | 蒽 (anthracene) | 3 | 14 | c1ccc2cc3ccccc3cc2c1 |
| 2 | 菲 (phenanthrene) | 3 | 14 | c1ccc2c(c1)ccc3ccccc32 |
| 3 | 菲变体 | 3 | 14 | c1ccc2c(c1)ccc3ccccc23 |
| 4 | 芘 (pyrene) | 4 | 16 | c1cc2ccc3cccc4ccc(c1)c2c34 |
| 5 | 三亚苯 (triphenylene) | 4 | 16 | c1ccc2c(c1)cc3cccc4ccc2c34 |
| 6 | 荧蒽 (fluoranthene) | 3 | 14 | c1ccc2c(c1)ccc3c2cccc3 |
| 7 | 吖啶 (acridine) | 3 | 14 | c1ccc2c(c1)nc3ccccc3c2 |
| 8 | 苯并[h]喹啉 | 3 | 14 | c1ccc2c(c1)cc3cccnc3c2 |
| 9 | 苯并[f]喹啉 | 3 | 14 | c1ccc2c(c1)ccc3cccnc23 |
| 10 | 苝 (perylene) | 5 | 20 | c1cc2cccc3ccc4cccc5ccc1c2c3c54 |
| 11 | 苯并[a]芘 | 4 | 18 | c1ccc2c(c1)cc3c4ccccc4ccc3c2 |

### 文献依据

- **Deep Lead Optimization (JACS 2024, §3.1-3.3)**: Scaffold diversity is the key to successful lead optimization. BM scaffold analysis shows 96% of drugs contain ring structures.
- **Round 18-19 实验验证**: Vina scoring function 偏向疏水芳香堆积。H025 已证明极性取代基降低 Vina 评分，反之大π表面积应提升评分。
- **MOOSE-Chem (Yang et al., 2025)**: "Diverse initial population is essential for evolutionary search to avoid premature convergence."

### 验证

- ✅ `python3 -m py_compile src/generator.py` 编译通过
- ✅ 所有81个骨架 SMILES 验证为有效分子
- ✅ 导入测试通过

### 预期效果

| 指标 | 基线 (H027) | 预期 |
|------|-------------|------|
| Best BE | -9.332 | -9.5 ~ -10.0 |
| Avg BE | -8.810 | -9.0 ~ -9.2 |
| Trivial ratio | 0/10 | 0/10 (无变化) |
| QED | ~0.45 | ~0.40 (略降，大芳香体系 MW 增加) |
