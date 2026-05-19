# 文献分析报告 — Round 10

> 基于 `papers/` 中三篇参考论文的深度分析，结合 9 轮迭代后的当前项目状态。

---

## 一、论文核心方法摘要

### 1. MOOSE-Chem (Yang et al., 2025) — 进化假设生成
- **核心机制**: 进化算法导航组合化学空间，评估函数引导选择
- **关键洞察**: "Diverse initial population is essential for evolutionary search to avoid premature convergence"
- **算子组合**: 变异 (mutation) + 重组 (crossover/recombination) + 骨架替换 (scaffold hopping)
- **本项目映射**: H001-H002 进化框架 + H010 scaffold hopping + H013 crossover → 框架已具备，但重组算子利用率低

### 2. MolLEO (Wang et al., 2024b) — LLM 作为进化算子
- **核心机制**: LLM 同时充当变异和重组算子，文本指令驱动
- **关键洞察**: 多种进化算子组合使用比单一算子更有效
- **本项目映射**: 当前 generator 只有 mutate/combine/random 三种策略，combine 策略实现有缺陷

### 3. ChemCrow (Bran et al., 2024) — 18 工具集成
- **核心启示**: 工具/规则库的丰富度直接决定 Agent 的能力边界
- **本项目映射**: synthesis_v2.py 已有 35+ 规则，trivial route 降至 0/10（此方向接近饱和）

### 4. Deep Lead Optimization (JACS 2024) — BRICS 碎片化
- **BRICS 分解**: 16 种可断裂键类型是分子碎片化的标准工具
- **四个核心子任务**: Scaffold Hopping ✓ (H010), Linker Design (部分), Side-chain Decoration (部分), Fragment Replacement ✗
- **关键启示**: Fragment replacement 是我们尚未实现的第四个核心子任务

---

## 二、技术要点与项目映射（更新版）

| 论文方法 | 技术要点 | 本项目映射 | 状态 |
|---------|---------|-----------|------|
| MOOSE-Chem 进化 | 进化算法 + 多样性 | generator.py | 已实现 |
| MOOSE-Chem 重组 | 分子片段重组算子 | generator.py combine/crossover | **combine 策略有 bug** |
| MolLEO 多算子 | 变异+重组+骨架替换 | generator.py | crossover 仅 15% |
| JACS Fragment Repl. | BRICS 片段替换算子 | — | **未实现** |
| JACS Side-chain Decor. | 骨架约束的侧链生长 | generator.py scaffold hopping | 部分实现 |
| LARC 逆合成 | 规则覆盖率 | synthesis_v2.py 35+规则 | 接近饱和 |
| JACS SA 阈值 | SA 2-5 为优质先导 | evaluator.py SA=6.0 | 已实现 |

---

## 三、按优先级排序的改进机会

### 高优先级（影响大 + 实现成本低）

| # | 假设ID | 改进方向 | 论文支撑 | 预期影响 |
|---|--------|---------|---------|---------|
| 1 | **H018** | 修复并增强 BRICS 片段重组策略 (combine) | MOOSE-Chem: 重组算子; JACS: BRICS 分解+重组; MolLEO: 多算子组合 | 生成全新化学型，突破当前骨架库限制 |
| 2 | H019 | 提升 Crossover 算子使用率（15%→30%） | MOOSE-Chem: 重组是进化核心算子 | 增加后代多样性 |

### 中优先级

| # | 假设ID | 改进方向 | 论文支撑 |
|---|--------|---------|---------|
| 3 | H020 | Fragment replacement 算子 (JACS 第四子任务) | JACS 2024: 片段替换是先导化合物优化的核心操作 |
| 4 | H021 | SNAr/卤素交换逆反应规则 | LARC: 规则覆盖率 |

---

## 四、本轮迭代方向

基于以上分析，本轮（Round 10）聚焦 **H018：修复 BRICS 片段重组 "combine" 策略**。

理由：
1. **当前 bug**：`strategy="combine"` 直接字符串拼接 SMILES（如 `c1ccccc1C(=O)Nc1ccccc1`），化学上无效
2. **高影响**：正确的 BRICS 重组可以从已有分子创造全新化学型，突破当前 SCAFFOLDS 库的局限
3. **小改动**：只需重写 combine 生成逻辑（~60 行），不涉及其他模块
4. **文献支撑强**：MOOSE-Chem + MolLEO + JACS 三方支撑
5. **时机恰当**：逆合成覆盖率已达 100%（trivial=0/10），结合能已接近天花板，分子多样性是当前最大增长空间
