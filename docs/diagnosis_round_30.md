# 诊断报告 — Round 30

> 时间: 2026-05-21
> 状态: 非首次运行，跳过完整文献解析，基于 knowledge_base.md 诊断

---

## 1. 当前基线 (H029 VERIFIED)

| Metric | Value |
|--------|-------|
| Best BE | **-10.099** kcal/mol |
| Avg BE | **-9.412** kcal/mol |
| Trivial ratio | **2/10** ⚠️ |
| Dominant chemistry | 100% 4-ring polycyclic aromatics (phenanthrene/pyrene derivatives) |
| LogP range | 3.39 - 4.91 |

## 2. 外部知识搜索

搜索关键词:
1. `Vina docking scoring function hydrophobic bias correction` → Google Scholar 返回多篇论文证实 Vina 存在系统性疏水偏倚
2. `TYK2 JH2 inhibitor binding affinity optimization` → 无明显新突破
3. `composite scoring logP penalty multi-objective docking` → 大量文献支持多目标评分修正

**关键发现**: Vina 疏水偏倚已被多篇文献记录 (JCIM 2018, 2020; J. Comput. Chem. 2010)。Tanchuk et al. (2016) 提出 hybrid scoring function 修正此偏差。

## 3. 瓶颈诊断

### 瓶颈 #1: Vina 疏水偏倚导致种群单一化

**根因**: Vina 评分函数中疏水项 (hydrophobic term) 对多环芳烃过度奖励。H029 多构象对接放大了此效应——更多构象采样让扁平分子找到更优姿态。

**证据**:
- Top-10 全部是 4 环多环芳烃，LogP 均值 4.25
- 无任何含 N/O 杂环铰链结合分子进入 Top-10
- 2 个 trivial 分子因异常稠环结构（环丙烷稠合、异噁唑稠合）无法被现有规则匹配

**影响**: 
- 化学多样性极低（单一化学型）
- Trivial ratio 退化（异常稠环无法合成）
- 所有分子都是扁平疏水型——与激酶铰链结合（需极性氢键）矛盾

### 瓶颈 #2: 复合评分缺少亲脂性惩罚

**根因**: 当前复合评分 `0.8×BE_norm + 0.2×route_quality` 对 LogP 无感知。
高 LogP 多环芳烃虽然 route_quality=0（trivial）但 BE 足够好，仍能胜出。

**证据**: 如果 BE 差距 > 25% 归一化范围，route_quality 的 0.2 权重无法扭转排名。

## 4. 改进假设

---

### 假设ID: H030
**瓶颈**: Vina 疏水偏倚 + 复合评分无 LogP 感知 → Top-10 被多环芳烃垄断
**文献支撑**: 
- Tanchuk et al. (Chem. Biol. Drug Des., 2016): Hybrid scoring 修正 Vina 疏水偏倚
- Gaillard (JCIM, 2018): Vina CASF-2013 评估确认疏水项过拟合
- Jacobsson & Karlén (J. Chem. Inf. Model., 2006): "Ligand bias of scoring functions" — 需要分子性质归一化

---

### ───────────────── 推理链 ─────────────────

| 步骤 | 内容 | 置信度 | 推理方式 | 依据来源 |
|------|------|--------|----------|----------|
| S1 | Vina 疏水项系统性过奖励高 LogP 多环芳烃 | 高 | 文献 | JCIM 2018, Tanchuk 2016 |
| S2 | 当前 Top-10 全部为 4-ring 多环芳烃 (LogP 3.4-4.9) | 高 | 实验 | output/result.csv |
| S3 | 不加 LogP 惩罚 → 多环芳烃始终主导排名 | 高 | 演绎 | 从 S1+S2 |
| S4 | 加 LogP 惩罚后排名会向低 LogP 分子倾斜 | 中 | 演绎 | 需验证 low-LogP 分子 BE 是否合理 |
| S5 | LogP 惩罚不会损害最佳 BE（对 extreme 分子保留） | 中 | 演绎 | 惩罚权重仅 10%，不影响极端高分分子 |

**综合置信度**: 中

S4/S5 为"中"——不确定低 LogP 分子的 BE 是否仍有竞争力。如果低 LogP 分子 BE 显著更差（如 >2 kcal/mol 差距），LogP 惩罚可能无法扭转。

---

### ───────────────── 验证标准 ─────────────────

**Q1: 如果核心指标提升 < 5%，是否仍保留？**
答: 是。此改进目标是化学多样性 + trivial ratio，而非结合能。

**Q2: 如果指标下降，最可能的原因是什么？**
答: 低 LogP 分子 BE 显著差于高 LogP 分子 → LogP 惩罚把非最优分子推上排名，Avg BE 可能退化。

**Q3: 本假设的最低可接受结果是什么？**
答: Trivial ratio 降至 1/10 以下，且化学型多样性增加（不再 100% 多环芳烃）。Avg BE 允许回落 0.5 kcal/mol。

---

### ───────────────── 改进方案 ─────────────────

**改动文件**: `tools/pipeline.py`
**改动内容**: 在复合评分中加入 LogP 惩罚项:
- `composite = 0.7×BE_norm + 0.2×route_quality + 0.1×(1 - logp_penalty)`
- `logp_penalty = max(0, (logP - 4.0) / 3.0)` — LogP ≤ 4 无惩罚，LogP=5 惩罚 33%, LogP≥7 惩罚 100%
- 需从 scored_candidates 获取 LogP 值

**验证指标**: pipeline 跑 50 分子，对比改前改后的:
- 最佳结合能 / Top-10 平均结合能
- Trivial route 比例
- 化学型分布（unique Murcko scaffolds 数量）
- 平均 LogP

