# 文献解析报告 — Round 1

## 目标蛋白
- **TYK2 激酶** (pdb|5C01|A, NON-RECEPTOR TYROSINE-PROTEIN KINASE TYK2)
- 257 aa, 对接坐标验证通过 (偏移 3.78 Å)

## 核心论文方法摘要

### 1. MOOSE-Chem (Yang et al., 2025) — 进化算法驱动的假设生成
**技术要点**：
- 将分子生成建模为进化优化问题：初始种群 → 变异/重组 → 适应度评估 → 选择
- **多目标优化**：同时优化 docking score 和 synthetic accessibility
- "Diverse initial population is essential for evolutionary search to avoid premature convergence"
- 使用 LLM 作为变异和重组算子

**映射到本项目**：
- H002 已实现 docking guidance（单目标 BE）
- H011 已实现 MMD 多样性保持
- **缺失**：SA score 仅用于过滤，未纳入适应度共优化

### 2. LARC (Baker et al., 2025) — Agent-as-a-Judge 逆合成
**技术要点**：
- 规则覆盖率是逆合成质量的关键决定因素
- 提出 Agent 自我评审路线可行性的机制
- 多步路线比单步路线更可靠（反映了真实合成复杂度）

**映射到本项目**：
- H003 已实现递归多步逆合成 (max_depth=3)
- 35+ 条逆合成规则
- **缺失**：路线质量评分机制 — 当前路线生成后无质量评估

### 3. Deep Lead Optimization (JACS, 2024) — 先导化合物优化四子任务
**技术要点**：
- **Scaffold Hopping**：替换核心骨架保留侧链
- **Linker Design**：优化片段间连接子
- **Fragment Replacement**：替换功能基团增强结合
- **Side-chain Decoration**：修饰侧链优化性质
- SA score 在 2-5 为优质先导化合物

**映射到本项目**：
- H010 已实现 Scaffold Hopping (20% 概率)
- **缺失**：Fragment Replacement 算子 — 可进一步提升分子多样性

### 4. ChemCrow (Bran et al., 2024) — 18工具集成
**技术要点**：
- 工具丰富度直接决定 Agent 探索的化学空间边界
- 集成 18 个专家化学工具
- 每个工具有明确的输入/输出接口

**映射到本项目**：
- 已有 6 个化学工具（生成/对接/合成/评估/pipeline/report）
- 逆合成规则库可进一步扩充

### 5. Coscientist (Boiko et al., 2023) — 自主实验
**技术要点**：
- 多 LLM 协作：Planner + Web Searcher + Code Execution + Automation
- "performing experiments multiple times" → 共识机制
- 基于实验结果的迭代反思

**映射到本项目**：
- H009 已实现共识对接（3次独立对接取中位数）
- ✅ 核心模式已采用

## 改进机会排序 (影响 × 易实现)

| 优先级 | 方向 | 影响 | 易实现 | 文献 |
|--------|------|------|--------|------|
| ⭐1 | 逆合成路线质量评分 + 规则优化 | 高 | 高 | LARC, ChemCrow |
| 2 | SA 共优化 fitness | 中 | 中 | MOOSE-Chem |
| 3 | Fragment Replacement 算子 | 中 | 中 | DLO JACS 2024 |
| 4 | 更多逆合成规则 | 低 | 高 | ChemCrow |

## 当前基线
- 最佳 BE: -9.941 kcal/mol, 平均 -8.961
- Trivial route: 1/10 (10%)
- 已实现 H001-H011
