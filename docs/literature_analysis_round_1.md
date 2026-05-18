# 文献解析报告 — Round 1

## 1. 论文核心方法摘要

### 1.1 ChemCrow (Bran et al., 2024) — 18 工具集成化学 Agent
- **技术要点**: LLM + 18 个专家设计化学工具（分子搜索、性质预测、反应预测等）
- **架构**: Toolbox-Based Tool Use — Agent 根据任务自动选择合适的工具链
- **对本项目的启示**: 工具/库的丰富度直接决定 Agent 探索的化学空间边界。本项目 SCAFFOLDS 库（55个骨架）和 RETRO_RULES（35条规则）的丰富度是核心竞争力

### 1.2 ChemAgents (Song et al., 2025) — 分层多 Agent 化学研究
- **技术要点**: Manager Agent 协调 5 个 Specialist Agent（文献阅读、机器人操作、数据分析等）
- **架构**: Hierarchical Tool Use in Multi-Agent Systems — 高层规划 + 专业执行
- **对本项目的启示**: 分而治之。当前 pipeline 是单线程的 "生成→对接→合成"，可考虑各阶段独立优化

### 1.3 MOOSE-Chem (Yang et al., 2025) — 进化算法导航化学假设空间
- **技术要点**: 使用进化算法（变异+重组+适应度选择）自动生成并验证化学假设
- **关键声明**: "Diverse initial population is essential for evolutionary search to avoid premature convergence"
- **对本项目的启示**: **核心洞察** — 当前 `generate_with_docking_guidance` 仅按结合能排序选择种子，缺少多样性保持机制，可能导致过早收敛于局部最优

### 1.4 MolLEO (Wang et al., 2024b) — LLM 作为变异和重组算子
- **技术要点**: LLM 直接充当进化算法中的变异和重组算子，基于文本指令的多目标优化
- **对本项目的启示**: 多目标优化（结合能 + QED + 多样性）优于单目标仅优化结合能

### 1.5 Deep Lead Optimization (Zhang et al., JACS 2024) — 先导化合物优化综述
- **技术要点**: 
  - 四类核心子任务：Scaffold Hopping, Linker Design, Side-Chain Decoration, Fragment Replacement
  - "Scaffold hopping 替换核心骨架同时保留有利取代基"
  - 优质先导化合物 SA 在 2-5 之间，QED 是药物相似性的综合指标
- **对本项目的启示**: H010 已实现 Scaffold Hopping；H001 已收紧 SA 阈值。可进一步在 Linker Design 方面增强

### 1.6 LARC (Baker et al., 2025) — Agent-as-a-Judge 逆合成
- **技术要点**: 规则覆盖率是逆合成质量的关键决定因素，多步路线评审提升规划质量
- **对本项目的启示**: H003 已实现递归多步逆合成。当前 trivial route 比例为 0%，说明合成规则覆盖面已足够

### 1.7 Coscientist (Boiko et al., 2023) — 自主化学实验
- **技术要点**: 
  - Multi-LLM 系统：Planner + Web Searcher + Code Execution + Automation
  - 核心模式："performing experiments multiple times" — 多次独立实验取共识
  - 迭代反思：根据实验结果修正代码
- **对本项目的启示**: H009 已实现 Consensus Docking；迭代反思已通过进化 pipeline 实现

## 2. 方法 → 代码映射

| 论文方法 | 当前实现 | 差距与改进机会 |
|---------|---------|--------------|
| MOOSE-Chem 多样性保持 | 选择仅按结合能排序 | **H011**: 引入多样性保持选择（Tanimoto 距离奖励） |
| MolLEO 多目标优化 | 单目标（结合能）排序 | **H012 候选**: 复合评分 = w1×BE + w2×QED + w3×diversity |
| ChemCrow 工具丰富度 | 55 scaffolds + 35 rules | 可继续扩充但收益递减 |
| Deep Lead Opt. Linker Design | 仅 C/O/N 三种连接子 | **H013 候选**: 扩充 LINKERS 库（酰胺、磺酰胺等） |
| ChemReasoner DFT 验证 | 仅 Vina docking | 受限于计算资源，暂不可行 |

## 3. 改进机会优先级排序

按「影响 × 实现容易度」排序：

| 优先级 | 假说ID | 方向 | 影响 | 易实现 | 得分 |
|--------|--------|------|------|--------|------|
| ⭐1 | H011 | 多样性保持选择 | 高 | 高 | ⭐⭐⭐ |
| 2 | H012 | 多目标复合评分 | 高 | 中 | ⭐⭐ |
| 3 | H013 | 扩充 LINKERS 库 | 中 | 高 | ⭐⭐ |
| 4 | H014 | 扩充 RETRO_RULES | 低 | 高 | ⭐ |

**结论**: 优先实现 H011（多样性保持选择），因为：
- MOOSE-Chem 明确警告过早收敛是进化算法的核心风险
- 修改量小（仅修改选择函数），风险低
- 可量化验证（多样性指标 + 结合能对比）
