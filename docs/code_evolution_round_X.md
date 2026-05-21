# 代码演进日志 — Round X (2026-05-21)

## 假设ID: H029 — 多构象对接增强

### 修改文件

1. **`src/docking.py`** — 核心修改
   - 新增函数 `smiles_to_pdbqt_multi_conformer(smiles, n_conformers=5)`：
     - 使用 RDKit `EmbedMultipleConfs` + ETKDGv3 生成多个 3D 构象
     - RMSD 剪枝阈值 0.5Å 保持构象多样性
     - MMFF94 优化每个构象
     - Meeko 转换为 PDBQT 文件列表
   - 修改 `dock_molecule()`：
     - 新增 `n_conformers` 参数（默认 1，保持向后兼容）
     - n_conformers > 1 时：对每个构象独立对接，取最优 BE
     - 返回新增 `n_conformers_tried` 和 `best_conformer_idx` 字段
   - 修改 `batch_dock()`：
     - 新增 `n_conformers` 参数，透传到 `dock_molecule`
   - 修改 `dock_molecule_consensus()`：
     - 新增 `n_conformers` 参数，透传到每次 `dock_molecule` 调用
   - 更新模块文档字符串：添加 H029 文献依据

2. **`tools/pipeline.py`** — 集成修改
   - 共识对接阶段启用 `n_conformers=3`（3 seeds × 3 conformers = 9 dockings/molecule）
   - 日志信息更新，标注 H009+H029 组合

### 文献依据

- Coscientist (Boiko et al., 2023): "performing experiments multiple times" 消除随机偏差
- GNINA Benchmarking (Molecules, 2025): 构象采样质量直接影响对接精度
- 计算化学基本原则: 构象系综对接优于单构象

### 设计决策

- 仅在共识对接阶段使用多构象（top-20 分子），进化阶段保持单构象
- 理由：平衡计算成本与精度提升。进化阶段需要高通量，共识阶段需要高精度
- 保持完整向后兼容：默认 n_conformers=1，不影响任何现有调用

### 预期效果

- Best BE 提升 0.2~0.5 kcal/mol（通过发现更好的结合模式）
- 共识对接 std 可能增大（构象间差异 > seed 间差异）
- 时间成本：共识阶段增加 ~2x（3 构象 × 3 seeds vs 1 构象 × 3 seeds）
