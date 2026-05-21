---
name: molcraft-paper-analysis
description: 文献解析关注点。仅首次运行/陷入困境需要回原文查找时读取。包含参考论文中值得关注的具体章节和概念。
---

# 文献解析指南

首次运行或陷入困境后强制恢复时使用。

## 综述论文关注点 (papers/autonomous_agents_survey.md)

### Chemistry Agent 部分
- **ChemCrow**: 18个工具集成的化学 Agent
- **ChemAgents**: 分层多 Agent（Manager+Specialist）
- **ChemReasoner**: LLM + DFT 假设验证
- **LARC**: Agent-as-a-Judge 逆合成
- **MOOSE-Chem**: 自动假设生成
- **FROGENT**: 端到端药物设计

### Multi-Agent 部分
- **TAIS**: 模拟研究团队
- **Agent Laboratory**: 文献 → 实验 → 论文

### Self-Code Evolution
- **AI Scientist**: 自主代码生成与迭代

## 解析产出

1. `docs/literature_analysis_round_X.md` — 完整分析报告
2. `docs/knowledge_base.md` — 策略库（每条: 标题 + 来源 + 技术要点 + 适用场景 + 已尝试标记）

## knowledge_base.md 格式

每条策略包含: 标题, 来源论文, 技术要点, 适用场景, 尝试状态(VERIFIED/REJECTED/PENDING)
Agent 每轮验证后更新标记。后续轮次直接读此文件，不必重读论文全文。