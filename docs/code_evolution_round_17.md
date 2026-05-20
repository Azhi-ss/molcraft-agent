# 代码演进 — Round 17 (H021)

## 假设 H021: TYK2 铰链结合骨架偏置

### 修改文件
- `src/generator.py`

### 修改内容

1. **新增 `KINASE_HINGE_SCAFFOLDS` 列表** (行 138-165):
   - 8 个激酶铰链结合特异性骨架:
     - 7-azaindole (吡咯并[2,3-b]吡啶)
     - 嘌呤 (两种变体)
     - 吡唑并[1,5-a]嘧啶
     - 咪唑并[1,2-a]吡啶
     - 喹唑啉
     - 吲哚
     - 嘧啶

2. **修改 `generate_molecules` — mutate 策略** (行 849-863):
   - 当 `scaffold is None` 时，每次种子选择以 40% 概率从 KINASE_HINGE_SCAFFOLDS 采样
   - 60% 概率从完整 SCAFFOLDS 采样（保持多样性）
   - 有显式 scaffold 参数时不受影响（保留精确控制能力）

3. **修改 `generate_molecules` — combine 回退** (行 891-900):
   - BRICS 重组产量不足时，回退变异也带激酶偏置

4. **修改 `generate_molecules` — random 策略** (行 907-916):
   - 同样加入 40% 激酶骨架偏置

### 文献依据
- MOOSE-Chem (Yang et al., 2025): "Targeted initial population design improves evolutionary convergence"
- JACS 2024: 骨架多样性 + 靶点特异性偏置提升命中率和分子质量
- Coscientist (Boiko et al., 2023): 基于领域知识的定向探索优于随机搜索
- PDB 5C01: TYK2 Met978 铰链区是经典 hinge binder 靶点

### 预期效果
- 激酶铰链骨架分子比例: 0/10 → ≥2/10
- 化学空间扩展: 从联芳/THIQ 扩展到嘌呤/氮杂吲哚等新化学型
- 结合能: 预期持平或改善（新化学型可能发现更好的结合模式）
