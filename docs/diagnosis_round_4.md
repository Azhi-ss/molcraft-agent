# 诊断报告 — Round 4

## 外部知识搜索摘要

- Semantic Scholar API: 429 rate-limited
- ChemRxiv: CloudFlare blocked
- PubMed: 403 forbidden
- **备选来源**: 本地论文库 `papers/deep_lead_optimization_jacs.md`（JACS 2024）、`papers/coscientist.md`（Boiko et al., 2023）

从本地论文中提取的关键洞察：
1. **JACS 2024 §3.1 Scaffold Hopping**: "replacing core structures while retaining functional activity" — 骨架多样性是先导化合物优化的基石
2. **JACS 2024 §3.3 Side-Chain Decoration**: scaffold anchors binding pattern, side chains enhance effects — 优先扩展骨架库而非侧链
3. **Coscientist §Main**: "The Agent can use tools to browse the internet... use other LLMs for various tasks" — 多层LLM集成模式可用于后续改进

## 代码审查发现

### 1. 骨架库缺口（generator.py SCAFFOLDS, 第38-136行）

当前135+骨架覆盖：
- 单环芳烃(6)、多环芳烃(2: 萘)、含氮单杂环(4)、饱和杂环(8)、含氮稠环(12)、药物骨架(6)、酰胺/磺酰胺(6)、卤代(5)、饱和环-芳环稠合(6)、稠杂环(5)、桥环(4)、螺环(3)、扩展饱和杂环(3)

**缺失的大芳香体系：**
- 蒽 (anthracene) — 3环线性稠合
- 菲 (phenanthrene) — 3环角形稠合
- 芘 (pyrene) — 4环稠合
- 荧蒽 (fluoranthene) — 4环
- 苯并萘/四芬 — 更大π体系

**文献依据:** JACS 2024指出"scaffold diversity is key for lead optimization"且96%药物含环结构。Vina scoring function的疏水项奖励大π表面积。

### 2. 突变算子分布（generator.py _mutate_mol, 第197-235行）

五种算子：add(25%), replace(25%), remove(10%), linker(20%), scaffold_hop(20%)。缺少环扩展/收缩算子（JACS 2024 scaffold hopping子类型），但实现复杂度高，本轮不优先。

### 3. 逆合成 / 评估 / 对接

- synthesis_v2.py: 53+规则，trivial 0/10 — **已饱和，非瓶颈**
- evaluator.py: 过滤器适当宽松(H027已调优) — **稳定**
- docking.py: Vina共识对接 — **稳定**
- pipeline.py: 2代进化 + 复合排序 — **稳定**

## 选定的改进假设

**H026 — Expanded Large Aromatic Scaffold Library**

| 维度 | 内容 |
|------|------|
| 瓶颈 | Scaffold库缺少3+环的大芳香体系（仅萘2变体） |
| 文献支撑 | JACS 2024: scaffold diversity + Vina hydrophobic bias |
| 改动范围 | generator.py SCAFFOLDS 列表 → 仅添加新条目 |
| 风险 | 低（纯增量，无现有逻辑修改） |
| 预期效果 | Best BE 从 -9.332 提升至 -9.5~-10.0, Avg BE 提升 3-5% |
| 验证指标 | Best BE, Avg BE, trivial ratio, QED |
