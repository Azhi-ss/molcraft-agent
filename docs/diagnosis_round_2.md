# 诊断报告 — Round 2 (H035)

## 1. 禁重检查

REJECTED: H019, H025, H026 → 均已排除。H035 方向不与任何 REJECTED 重叠。

## 2. 外部知识搜索

本轮聚焦代码审查，未进行额外外部搜索（Round 1 已搜索 LKM + arXiv）。

## 3. 源码审查 — 关键发现

### 多构象对接 (H029) 已实现但未启用

`src/docking.py`:
- `smiles_to_pdbqt_multi_conformer()` (line 70): 支持生成 n 个 ETKDGv3 构象 ✅
- `dock_molecule(n_conformers=1)` (line 140): 默认 n_conformers=1 ❌
- `batch_dock(n_conformers=1)` (line 310): 默认 n_conformers=1 ❌

`tools/pipeline.py`:
- 三处 `batch_dock()` 调用均未传 `n_conformers` → 使用默认 1

### GNINA 基准文献结论

构象采样质量直接影响对接精度（Molecules, 2025）:
- 单构象对接可能错过最佳结合模式
- 多构象系综对接（3-5 构象）中位数提升约 0.3-0.6 kcal/mol
- 成本: 线性增加（n×对接时间）

## 4. 瓶颈诊断

**瓶颈**: 当前每次对接仅使用单一 3D 构象（ETKDGv3 单次嵌入），可能漏掉更优的结合构象。这解释了为什么 BE 长期稳定在 -9.6~-9.9 区间无法突破——对接算法已找到当前构象的局部最优，但更优的全局构象未被采样。

## 5. 假设提出

### H035: Enable Multi-Conformer Docking (n_conformers=3)

**来源**: GNINA Benchmarking (Molecules, 2025) + 已实现的 H029 多构象对接基础设施

**方案**:
1. 修改 `dock_molecule()` 默认 `n_conformers=1→3`
2. 修改 `batch_dock()` 默认 `n_conformers=1→3`  
3. 运行 pipeline 验证 BE 改善

**预期效果**: 
- Best BE 改善 0.3-0.6 kcal/mol（文献基准）
- Avg BE 相应改善
- 对接时间增至 3×（约 15-30 min → 45-90 min）

**验证标准**:
- Best BE 低于 -9.9 kcal/mol
- Avg BE 不低于 -8.5 kcal/mol
- Trivial 保持 0/10
- 对接成功率不下降
