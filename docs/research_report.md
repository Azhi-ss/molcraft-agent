# MolCraft Agent — 科研报告

> TYK2 (5C01) 抑制剂设计 | 自主科研迭代 | 5轮实验

---

## 1. 文献解析关键发现

从 `papers/` 目录中的参考论文中提取的核心洞察：

### Deep Lead Optimization (JACS 2024)
- **四大先导优化子任务**: Scaffold Hopping, Side-Chain Decoration, Linker Design, Fragment Replacement
- **BRICS 分解**: 16种可断裂键是片段重组的工业标准
- **关键指标**: SA score 2-5 为优质先导化合物范围，MW 150-500 为类药区间

### Coscientist (Boiko et al., 2023)
- **Multi-LLM Agent 架构**: Planner + Web Searcher + Code Executor + Automation
- **迭代反思**: "基于实验结果的迭代反思"是 Agent 持续优化的核心
- **共识验证**: "performing experiments multiple times" 消除随机噪声

### MOOSE-Chem (Yang et al., 2025)（通过综述转引）
- **进化算法**: 生成→对接→选择→变异循环是分子优化的最优框架
- **多样性保持**: MMD 贪心选择防止过早收敛
- **Crossover**: 片段重组（recombination）是核心进化算子

---

## 2. 诊断瓶颈与假设

### 已验证假设

| ID | 描述 | 结果 |
|----|------|------|
| H001 | SA score 阈值 8.0→6.0 | ✅ VERIFIED |
| H002 | Docking-guided generation | ✅ VERIFIED (+0.4~0.8 kcal/mol) |
| H003 | 递归多步逆合成 (max_depth=3) | ✅ VERIFIED |
| H009 | 共识对接 (3次取中位数) | ✅ VERIFIED |
| H010 | Scaffold Hopping 变异算子 | ✅ VERIFIED |
| H011 | MMD 多样性保持选择 | ✅ VERIFIED (+19.3% BE) |
| H012 | Route Quality Scoring (Agent-as-a-Judge) | ✅ VERIFIED |
| H013 | Crossover 重组算子 | ✅ VERIFIED |
| H014 | Mass Balance 化学计量验证 | ✅ VERIFIED |
| H015-H020 | 逆合成规则扩容 (53+规则) | ✅ VERIFIED (trivial 0/10) |
| H022 | 桥环/螺环骨架库扩展 | ✅ VERIFIED |
| H027 | 宽松分子过滤器 (max_rings 7→9) | ✅ VERIFIED (+2% BE) |
| H028 | Suzuki SMARTS 特异性修复 | ✅ VERIFIED |

### 已拒绝假设

| ID | 描述 | 失败原因 |
|----|------|----------|
| H019 | Suzuki Byproduct Atom Balance | Trivial ratio 退化 |
| H025 | 极性药效团取代基扩展 | Vina 惩罚极性基团 (-2% BE) |
| H026 | 大芳香环骨架扩展 (蒽/菲/芘) | 骨架未被选用，MW 超限 |

---

## 3. 代码演进修改摘要

### Round 4 (H026 — REJECTED)
- **generator.py**: 新增12个大芳香环骨架（蒽、菲、芘、苝、吖啶、苯并喹啉、苯并芘）
- **影响**: 无正面效果，被 MW=500 过滤拦截

### Round 5 (H028 — VERIFIED)
- **synthesis_v2.py**: Suzuki SMARTS `[c;R][c;R]` → `[c;R]!@[c;R]`
- **问题**: 原 SMARTS 错误匹配稠环体系内 C-C 键（如 indole 有13个假匹配）
- **修复**: `!@` 限定仅匹配非环键（即联芳基 C-C 键）
- **效果**: Trivial route 回归 1/10→0/10 修复

---

## 4. 实验验证结果

### 最终配置
```
n_generate=50, n_top=10, strategy=mutate, n_generations=2, docking_guidance=True
```

### 最终性能

| 指标 | 值 |
|------|-----|
| 最佳结合能 | **-9.19 kcal/mol** |
| 平均结合能 | **-8.42 kcal/mol** |
| Trivial route 比例 | **0/10** |
| 主导化学型 | Suzuki 联芳基 + 喹唑啉/咪唑并吡啶 |

### Top 3 分子

| # | SMILES | BE (kcal/mol) | 合成路线 |
|---|--------|---------------|----------|
| 1 | O=Cc1cc2cccc(C(=O)O)c2cc1-c1ccccc1 | -9.19 | 3步: 氧化→溴化→Suzuki |
| 2 | Cc1ccc(-c2cncc3ccccc23)cc1 | -9.11 | 2步: Pictet-Spengler→Suzuki |
| 3 | FCc1ncc2cccc(-c3ccccc3)c2n1 | -8.73 | 2步: 喹唑啉合成→Suzuki |

---

## 5. 迭代过程科学洞察

1. **Vina 评分偏向确认**: 疏水芳香堆积是 TYK2 ATP 口袋的主要驱动力，极性药效团反而降低评分（H025 验证），大芳香体系也无效（H026 验证）

2. **SMARTS 特异性至关重要**: 过宽的 SMARTS 模式会导致假匹配，阻断特定规则。`[c;R]!@[c;R]` 优于 `[c;R][c;R]` 是防止稠环假匹配的关键（H028 发现）

3. **分子过滤器需要精确校准**: H027 的宽松过滤器 (+2% BE) 和 H026 的失败表明，过滤器的"甜区"在 MW 250-400 区间，过大或过小都不利

4. **逆合成规则覆盖率已达饱和**: 53+规则实现 trivial 0/10，进一步扩充的边际收益递减

5. **进化框架成熟**: 2代进化+MMD多样性+共识对接+路线质量评分已形成稳健的端到端流程
