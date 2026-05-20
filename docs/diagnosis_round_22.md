# 诊断报告 — Round 22 (H026)

## 当前基线

| Metric | Value |
|--------|-------|
| Best BE | -9.332 |
| Avg BE | -8.810 |
| Trivial ratio | 0/10 |
| Dominant chemistry | Imidazopyridine + Suzuki biaryl |

## 瓶颈分析

### 瓶颈 1: 大芳香骨架库不足

**现状**: SCAFFOLDS 列表 63 个骨架，最大芳香系统为萘（2环）。缺少数百种药物中常见的3-4环大芳香骨架（蒽、菲、芘等）。

**证据链**:
- H025 (极性取代基) 被拒绝：Vina 惩罚极性基团，偏好疏水平面芳香体系
- H027 (放宽过滤器 max_rings 7→9) 被验证：Avg BE +4.2%，证明更大分子有优势
- Deep Lead Optimization (JACS 2024): "the variety of parent skeletons should be notably broad" 骨架多样性决定化学空间覆盖
- Round 14 历史最佳 -10.359 来自扩展骨架库 (H022)，证明骨架多样性直接提升 BE

**根本原因**: Vina 评分函数将疏水堆积和 π-π 相互作用作为主要能量贡献项。更大芳香体系提供更多疏水表面积和 π 电子云重叠面积，Vina 会给予更优评分。

### 瓶颈 2: 合成路线化学类型单一

所有 top 分子均使用 Suzuki 偶联断开联芳基键。虽然路线有效（0/10 trivial），但缺乏化学多样性。不过这不是当前优先瓶颈——现有合成规则库已足够支撑当前分子类型。

## 搜索外部知识

⚠️ 网络不可用（PubMed、bioRxiv、ChemRxiv、Wikipedia 均 blocked）。使用本地论文 `papers/deep_lead_optimization_jacs.md` 作为知识来源。

**关键发现**:
- §3.2: BRICS 将分子分解为片段，16 种可断裂化学键
- §4.2 Scaffold Hopping: 骨架跃迁是核心先导化合物优化策略
- §2: "the prior distribution P_S(G\μ) should be highly expressive, or the variety of parent skeletons should be notably broad"
- 96% 的药物含环结构，56% 的分子量来自环组分

## 提出假设

### H026 — 扩展大芳香骨架库

```
假设ID: H026
瓶颈: SCAFFOLDS 列表缺少 ≥3 环的大芳香烃和杂芳烃骨架
文献支撑: Deep Lead Optimization (JACS 2024) — 骨架多样性决定化学空间探索效率;
          H022/H027 已分别验证扩展骨架和放宽过滤器的有效性；
          Vina 偏爱疏水芳香体系的实验证据 (H025 拒绝, H027 验证)

───────────────── 推理链 ─────────────────
步骤 | 内容                                  | 置信度 | 推理方式 | 依据来源
S1   | Vina 评分函数奖励疏水π-π堆积           | 高     | 实验     | H025 拒绝(极性→BE↓), H027 验证(大分子→BE↑)
S2   | 更大芳香系统提供更多疏水表面积          | 高     | 物化原理 | π电子离域面积∝环数
S3   | 添加蒽/菲/芘等骨架会增加高分分子概率     | 高     | 演绎     | 从 S1+S2
S4   | 现有合成规则能处理这些骨架的断键         | 中     | 演绎     | Suzuki断联芳基已覆盖，但稠环体系可能需要Diels-Alder或Friedel-Crafts规则

综合置信度 = 中（S4 的断键覆盖是不确定因素）
若验证失败，最可能原因是新骨架生成的分子的逆合成路线回归 trivial。

───────────────── 验证标准 ─────────────────
Q1 如果核心指标提升 < 5%，是否仍保留？
答：是 — 只要不退化即可保留（扩充骨架库本身就有长期价值）

Q2 如果指标下降，最可能的原因是什么？
答：大芳香分子 SA score 过高被过滤器拒绝，或逆合成路线回退 trivial

Q3 本假设的最低可接受结果是什么？
答：Best BE ≥ -9.15 (不退化)，Trivial ≤ 1/10

───────────────── 改进方案 ─────────────────
改动文件: src/generator.py
改动内容: SCAFFOLDS 列表新增 10-15 个大芳香骨架：
  - 蒽 (anthracene)
  - 菲 (phenanthrene)  
  - 芘 (pyrene)
  - 芴 (fluorene)
  - 荧蒽 (fluoranthene)
  - 吖啶 (acridine)
  - 吩嗪 (phenazine)
  - 吩噻嗪 (phenothiazine)
  - 咔唑 (carbazole)
  - 二苯并呋喃 (dibenzofuran)
  - 二苯并噻吩 (dibenzothiophene)
  - 苯并[a]蒽变体
验证指标: Best BE、Avg BE、Trivial ratio；至少1个分子含≥3环芳香体系
```
