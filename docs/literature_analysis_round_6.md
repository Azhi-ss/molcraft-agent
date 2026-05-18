# 文献解析报告 — Round 6

## 论文核心方法摘要

### 1. LARC: Agent-as-a-Judge 逆合成框架 (Baker et al., 2025)
- **技术要点**: 将 LLM 用作逆合成路线的"评审员"，对每条候选路线做多维度评估（化学计量守恒、反应可行性、原料可得性）
- **关键洞察**: 规则覆盖率是逆合成质量的决定因素，但同等重要的是对规则产出的验证——不匹配的路线应被拒绝而非标记为 trivial
- **映射到本项目**: synthesis_v2.py 已有 score_route_quality 评审函数（H012），但缺少前置的化学计量验证环节 — 规则匹配成功后应验证产物=原料之和

### 2. MOOSE-Chem: 进化算法导航组合空间 (Yang et al., 2025)
- **技术要点**: 进化算法（变异+重组+选择）作为分子生成的优化框架，多样性初始种群防过早收敛
- **关键洞察**: "Diverse initial population is essential for evolutionary search" — 已通过 H011 多样性选择实现
- **映射到本项目**: generator.py 已有成熟的进化生成框架，无需进一步优化

### 3. Deep Lead Optimization (JACS 2024)
- **技术要点**: 先导化合物优化四个子任务：Scaffold Hopping、Linker Design、Side-chain Decoration、Fragment Replacement
- **关键洞察**: 每个子任务对应特定的化学转化逻辑，而非随机变异
- **映射到本项目**: Scaffold Hopping (H010) 和 Side-chain decoration 已实现；Linker Design 目前是随机插入连接子，可优化

### 4. ChemCrow: 18工具集成化学Agent (Bran et al., 2024)
- **技术要点**: 工具/规则库的丰富度直接决定 Agent 能力边界
- **映射到本项目**: synthesis_v2.py 的 REACTION_RULES 已有 35+ 条规则

## 改进机会（按影响×易实现排序）

| 优先级 | 改进方向 | 影响 | 难度 | 文献依据 |
|--------|----------|------|------|----------|
| ⭐1 | 逆合成化学计量验证（H014）| 高 | 低 | LARC 2025 |
| 2 | 取代基保留的逆合成规则 | 中 | 中 | JACS 2024 |
| 3 | 稠环/桥环断键规则 | 中 | 高 | LARC 2025 |
