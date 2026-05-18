# 代码演进日志 — Round 2 (H013)

## 假设 H013: 分子 Crossover 重组算子

### 修改文件

#### `src/generator.py`

**新增函数 `_crossover_mol(mol1, mol2)`**：
- Murcko 骨架交换算法：取 parent2 骨架 + parent1 侧链 → 新分子
- 连接点自动检测：在核心上找可用碳原子，侧链上找未饱和原子
- 最多连接 3 个侧链片段，避免过度复杂
- 验证：SanitizeMol + SMILES 往返 + 与父代不同
- 成功率 ~30%（进化算法中可接受）

**修改 `generate_with_docking_guidance`**：
- 后续代（gen > 0）中，15% 概率执行 crossover
- 从种子池随机选 2 个不同分子作为双亲
- 85% 概率执行原有单亲变异

### 文献依据
- MOOSE-Chem (Yang et al., 2025): "Evolutionary operators include crossover between parent molecules, fragment swapping, and scaffold hopping"
- MolLEO (Wang et al., 2024b): LLM 驱动的重组操作提升化学空间探索效率
- Deep Lead Optimization (JACS 2024): Side-chain decoration + Scaffold Hopping 可组合产生新化学型

### 编译验证
- `src/generator.py`: ✅ 编译通过
- `_crossover_mol` 测试: ✅ 3/10 成功率，产生合理分子
