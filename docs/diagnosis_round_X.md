# 诊断报告 — Round X (2026-05-21)

## 当前基线

| Metric | Value |
|--------|-------|
| Best BE | -9.19 kcal/mol |
| Avg BE | -8.42 kcal/mol |
| Trivial ratio | 0/10 |
| 主导化学型 | Suzuki biaryl + quinazoline/imidazopyridine scaffolds |

## 外部知识搜索结果

Scholar 搜索获取的关键文献：
1. **GNINA vs Vina Benchmarking** (Molecules, 2025): GNINA 使用 CNN-based scoring，在 VS 场景中显著优于 Vina
2. **Vinardo Scoring Function**: AutoDock Vina 内含的替代评分函数，参数化不同
3. **GSScore** (Briefings in Bioinformatics, 2024): Graphormer-based shell-like scoring
4. **XDock** (2024): General docking method for protein-ligand interactions

关键洞察：Vina 评分函数已知偏向疏水芳香体系（与 H025/H027 发现一致），构象采样质量直接影响对接精度。

## 瓶颈诊断

### 瓶颈 1: 单构象对接限制 (优先级: HIGH, 影响: MEDIUM)
- **位置**: `src/docking.py` — `smiles_to_pdbqt()` 函数 Line 24-27
- **问题**: `EmbedMolecule(mol, randomSeed=42)` 仅生成单构象，MMFF 优化后直接对接
- **根因**: 分子可能有多个低能构象，单构象对接可能错过最优结合模式
- **文献支撑**: 
  - Coscientist (Boiko 2023): "performing experiments multiple times"
  - GNINA 论文 (2025): 强调构象采样对对接精度的影响
  - 药物设计基本原则: 构象系综对接优于单构象

### 瓶颈 2: Vina 评分系统性偏向 (优先级: MEDIUM)
- **位置**: `src/docking.py` Line 56 — `Vina(sf_name='vina', ...)`
- **问题**: Vina 默认评分函数对疏水芳香体系的过度奖励 (H025 已验证)
- **替代方案**: Vinardo 评分函数 (Vina RaDii Optimized) 可尝试

### 瓶颈 3: 逆合成规则覆盖率饱和 (优先级: LOW)
- **位置**: `src/synthesis_v2.py` — `RETRO_RULES`
- **问题**: 53+ 规则已实现 trivial 0/10，边际收益递减
- **判定**: 当前阶段无需扩展

## 提出假设

### H029 — 多构象对接增强 (PRIORITY: HIGHEST)

```
假设ID: H029
瓶颈: docking.py 仅生成单构象进行对接，可能错过最优结合模式
文献支撑: Coscientist (Boiko 2023) — "performing experiments multiple times";
          GNINA Benchmarking (Molecules 2025) — 构象采样质量影响对接精度

───────────────── 推理链 ─────────────────
步骤 | 内容                                    | 置信度 | 推理方式 | 依据来源
S1   | 单构象无法代表分子的完整构象空间           | 高     | 基本原理 | 计算化学
S2   | 多个构象独立对接可探索更多结合模式         | 高     | 演绎     | 从 S1
S3   | 从多构象结果中取最优 BE 可能优于单构象     | 高     | 演绎     | 从 S2
S4   | 5x 计算开销可接受 (每个分子多 20-30s)     | 中     | 估计     | 经验

综合置信度 = 中 (S4 为"中" — 时间成本需要实验验证)

───────────────── 验证标准 ─────────────────
Q1 如果核心指标提升 < 5%，是否仍保留？
答：是，理由：即使 BE 提升有限，多构象方法更符合科学最佳实践

Q2 如果指标下降，最可能的原因是什么？
答：多构象带来的额外自由度导致 Vina 评分噪声增加

Q3 本假设的最低可接受结果是什么？
答：Best BE 不低于基线 (-9.19) 的 95%（即 ≥ -8.73），Avg BE 不低于 -7.98

───────────────── 改进方案 ─────────────────
改动文件: src/docking.py
改动内容: 
  1. 新增函数 smiles_to_pdbqt_multi_conformer(smiles, n_conformers=5)
  2. dock_molecule() 中使用多构象：对每个构象对接取最优
  3. 保留原有单构象功能（向后兼容）
验证指标: pipeline 跑 50 个分子，2 代进化，对比改前改后 Best/Avg BE
```
