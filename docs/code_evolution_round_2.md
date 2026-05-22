# 代码演进报告 — Round 2 (H035)

## 修改文件

### `src/docking.py` — 多构象对接默认值更新

**修改内容**: 将三个关键函数的 `n_conformers` 默认值从 1 改为 3。

| 函数 | 旧默认值 | 新默认值 | 行号 |
|------|---------|---------|------|
| `dock_molecule()` | n_conformers=1 | **n_conformers=3** | 140 |
| `dock_molecule_consensus()` | n_conformers=1 | **n_conformers=3** | 257 |
| `batch_dock()` | n_conformers=1 | **n_conformers=3** | 310 |

**影响范围**:
- `dock_molecule()`: 所有单体对接（包括对接引导内循环）→ 每分子生成 3 个构象独立对接，取最优
- `dock_molecule_consensus()`: 共识对接 → 3 构象 × 3 种子 = 9 次独立对接，取中位数
- `batch_dock()`: 管线批量对接 → 所有分子使用多构象

**成本估算**:
- 每分子对接时间: 5-15s → 15-45s (3×)
- Pipeline 总耗时: ~10-20 min → ~30-60 min
- 内存/磁盘: 临时 PDBQT 文件增加 3×，自动清理不累积

**文献依据**:
- GNINA Benchmarking (Molecules, 2025): 构象采样质量直接影响对接精度
- Coscientist (Boiko et al., 2023): 重复实验消除随机偏差（多构象 = 构象空间的重复实验）
- H029 已实现基础设施，本修改只是启用

## 编译检查

✅ `python3 -m py_compile src/docking.py` 通过

## 风险评估

- **向后兼容**: 仅修改默认值，所有显式传参不受影响
- **性能**: 3× 对接时间，可接受（pipeline 原 12 min → 预计 36 min）
- **行为变化**: 对接结果可能因更优构象而改善，不会退化（取最优构象能量）
- **并发**: 不影响并发逻辑，仅增加每分子计算量
