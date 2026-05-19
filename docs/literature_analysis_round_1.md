# 文献分析报告 — 第1轮

## 靶点蛋白
- **名称**: TYK2 (Non-receptor Tyrosine-protein Kinase) — PDB: 5C01 Chain A
- **结构**: 257个氨基酸, 1条链, 激酶催化域
- **活性位点**: [19.7, 1.18, 24.76] (clft), 对接坐标偏移 3.78 Å（OK）

## 论文核心方法摘要

### 1. LARC — Agent-as-a-Judge 逆合成框架 (Baker et al., 2025)
- **技术要点**: 规则覆盖率是逆合成质量的核心决定因素；Agent 评审路线质量
- **映射到本项目**: 直接影响 `synthesis_v2.py` 的 `REACTION_RULES` 数量和 `score_route_quality()` 函数
- **已实现**: H012 路线质量评分、H014 质量守恒验证

### 2. ChemCrow — 18工具集成化学Agent (Bran et al., 2024)
- **技术要点**: 工具/知识库的丰富度直接决定Agent能力边界
- **映射到本项目**: REACTION_RULES 规则库大小、SCAFFOLDS 骨架库大小
- **已实现**: SCAFFOLDS 已扩充至 55 个骨架（H007）、RETRO_RULES 已扩充至 35+ 条（H012/H016）

### 3. MOOSE-Chem — 进化算法假设生成 (Yang et al., 2025)
- **技术要点**: 进化搜索导航组合空间；多样化初始种群防止过早收敛
- **映射到本项目**: `generate_with_docking_guidance()` 的进化循环、H011 多样性选择
- **已实现**: H002 对接引导生成、H011 MMD 多样性保持、H013 Crossover 重组

### 4. Coscientist — 多LLM自主实验 (Boiko et al., 2023)
- **技术要点**: Planner→Web Searcher→Code Execution 管道；多次实验取共识
- **映射到本项目**: H009 共识对接（3次独立对接取中位数）
- **已实现**: dock_molecule_consensus()

### 5. Deep Lead Optimization — 先导化合物优化四子任务 (JACS 2024)
- **技术要点**: Scaffold Hopping, Linker Design, Fragment Replacement, Side-chain Decoration
- **映射到本项目**: H010 Scaffold Hopping 算子、evaluator.py 的 SA score 阈值
- **已实现**: H001 SA 阈值收紧（8.0→6.0）、H010 骨架替换算子

## 改进机会列表（按影响×实现难度排序）

| 优先级 | 方向 | 影响 | 难度 | 状态 |
|--------|------|------|------|------|
| 1 | 扩充逆合成规则库（稀疏键类型） | 高 | 低-中 | ⬜ 本轮 |
| 2 | 引入 Enamine REAL 等商业库做起始原料验证 | 中 | 中 | ⬜ |
| 3 | 多步路线中引入保护基策略识别 | 中 | 高 | ⬜ |
| 4 | 对接盒子微调（向活性位点 cleft 偏移） | 低-中 | 低 | ⬜ |
| 5 | 引入 ChemAgents 分层架构 | 高 | 高 | ⬜ |

## 本轮重点
**扩充逆合成规则库**：当前 REACTION_RULES 约 35+ 条，但 lactone 开环、sp3C-sp3C 断键、环氧开环等常见转化仍缺失。本轮将基于 result.csv 中的 trivial route 分子分析断键缺口，设计 5-10 条新规则。
