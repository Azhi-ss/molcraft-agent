# 文献分析报告 — Round 8

> 基于 `papers/` 中三篇参考论文的系统分析，提取可落地的架构设计与方法论改进。

---

## 一、论文核心方法摘要

### 1. Coscientist (Boiko et al., 2023)
- **核心架构**: Planner + Web Searcher + Docs Searcher + Code Execution + Automation
- **关键创新**: 自纠错循环 — Agent 在代码出错后自动读取 traceback、修改代码、重新执行（Suzuki反应计算中SymPy缺失→basic Python→加print()的3轮自纠）
- **已验证能力**: 成功自主执行 Suzuki/Sonogashira 交叉偶联反应，从搜索文献→计算用量→编写协议→执行实验
- **安全机制**: 对已知危险化合物（海洛因、芥子气）进行拒止，但对未知化合物检测不足

### 2. Autonomous Agents Survey (Zhou et al., 2025)
- **ChemCrow (Bran et al., 2024)**: 18个专家设计化学工具集成，覆盖有机合成、药物发现、材料设计。**核心启示：工具/规则库的丰富度直接决定 Agent 的能力边界。**
- **ChemAgents (Song et al., 2025)**: 分层多Agent架构 — Manager协调多个Specialist（Literature Reader、Robot Operator等）。**核心启示：复杂工作流需要分工协作。**
- **MOOSE-Chem (Yang et al., 2025)**: 进化算法驱动假设生成，结合启发式搜索和背景知识验证。**核心启示：多样性初始种群对避免进化搜索过早收敛至关重要。**
- **LARC (Baker et al., 2025)**: Agent-as-a-Judge 逆合成，规则覆盖率是逆合成质量的关键决定因素。
- **AI Scientist (Yamada et al., 2025)**: 自主代码生成与迭代 — LLM生成实验代码、运行、分析结果、改进代码。**核心启示：闭环自进化是长期目标。**
- **TAIS (Liu et al., 2024)**: 模拟研究团队的多Agent协作框架。

### 3. Deep Lead Optimization (Zhang et al., JACS 2024)
- **四个核心子任务**: Scaffold Hopping、Linker Design、Side-chain Decoration、Fragment Replacement
- **关键数据**: SA score 2-5 是优质先导化合物的范围；BM Scaffold 是药用化学中最重要的骨架表示法
- **BRICS vs RECAP**: 16种可断裂键类型是分子碎片化的基础工具
- **方法论启示**: 约束子图生成统一视角 — de novo设计和lead optimization可以互补

---

## 二、技术要点与项目映射

| 论文方法 | 技术要点 | 本项目映射 | 当前状态 |
|---------|---------|-----------|---------|
| Coscientist 自纠错 | 读取traceback→修改→重试 | docking.py 已有异常处理 | 基础版已实现 |
| ChemCrow 18工具 | 工具丰富度决定能力边界 | synthesis_v2.py 35+规则 | 已实现，可继续扩充 |
| MOOSE-Chem 进化 | 进化算法+多样性保障 | generator.py scaffold hopping + MMD选择 | 已实现 H010/H011 |
| LARC 逆合成 | 规则覆盖率+多步规划 | synthesis_v2.py 递归规划 max_depth=3 | 已实现 H003 |
| JACS SA阈值 | SA 2-5为优质先导 | evaluator.py SA=6.0, rings≤7 | 已实现 H001 |
| JACS BRICS | 16种可断裂键 | BRICS fallback in synthesis_v2 | 已实现 |

---

## 三、按优先级排序的改进机会

### 高优先级（影响大 + 易实现）

| # | 假设ID | 改进方向 | 论文支撑 | 预期影响 |
|---|--------|---------|---------|---------|
| 1 | H016 | 添加 Diels-Alder 逆反应规则（稠环/桥环体系断键） | JACS 2024: 稠环体系需要专门断键策略; LARC: 规则覆盖率决定质量 | 降低BRICS fallback，提升路线合理性 |
| 2 | H017 | 添加 Michael/aldol/Claisen C-C键构建逆反应规则 | ChemCrow: 工具丰富度; LARC: 逆合成规则覆盖 | 覆盖更多非芳香C-C键断键 |
| 3 | H018 | 引入分子组合(combine)生成策略提升多样性 | MOOSE-Chem: 多样性种群; MolLEO: 多种进化算子 | 提升分子多样性，避免局部最优 |

### 中优先级

| # | 假设ID | 改进方向 | 论文支撑 |
|---|--------|---------|---------|
| 4 | H019 | 添加 SNAr 和杂环功能化逆反应规则 | LARC + JACS |
| 5 | H020 | 添加螺环特异性逆反应规则 | JACS: 螺环形成需针对性策略 |

### 低优先级（长期目标）

| # | 方向 | 说明 |
|---|------|------|
| 6 | 多Agent协作（ChemAgents风格） | 需要大规模架构改动 |
| 7 | 自纠错闭环（AI Scientist风格） | 需要基础设施改动 |

---

## 四、本轮迭代方向

基于以上分析，本轮（Round 8）聚焦 **H016：添加 Diels-Alder 环加成逆反应规则**。

理由：
1. 稠环体系在药物分子中极为常见（96%药物含环结构 — JACS 2024）
2. Diels-Alder 是构建六元环/桥环的最经典方法
3. 当前规则库没有覆盖此类断键
4. 实现只需 3-5 条 SMARTS 规则，修改量小

**下一轮（Round 9）**计划实施 H017（C-C键构建规则）或 H018（combine策略）。
