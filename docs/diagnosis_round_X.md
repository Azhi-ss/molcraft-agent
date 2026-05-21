# 诊断报告 — Round X (新会话)

## 日期: 2026-05-21

---

## 1. 禁重检查

REJECTED 假设:
- H019: Suzuki Byproduct Atom Balance — trivial ratio regression
- H025: Polar Pharmacophore Substituents — Vina penalizes polarity
- H026: Large Aromatic Scaffolds — significant BE regression

新假设方向均不与此三条重叠。✅

---

## 2. 当前基线 (result.csv 分析)

| 指标 | 值 |
|------|-----|
| 分子数 | 10 |
| Best BE | -9.972 |
| 实际 trivial route 数 | **2/10** (Row 6-7: OH→Cl 单原子替换) |
| 报告中 trivial ratio | 0/10 (错误) |

### Trivial 路线详情:
- Row 6: `O=C1c2ccccc2NNc2cccc(Cl)c21` → `...NNc2cccc(O)c21.Cl>>...NNc2cccc(Cl)c21.O`
  - 本质: 芳基 Cl→OH 单原子替换
  - `_is_trivial_route` 返回 False (reactant≠product)
  - `route_quality` = ~0.7 (未被识别为 trivial)

- Row 7: `CC1CCC2(CC1)Cc1cccc(Cl)c1C2` → `...Cc1cccc(O)c1C2.Cl>>...Cc1cccc(Cl)c1C2.O`
  - 本质: 同上

---

## 3. 瓶颈分析

### 瓶颈: Trivial 路线检测不完整

**根因**: `synthesis_v2.py` 第 407-409 行的卤素交换规则:
```python
("[c]Cl", "[c:1]Cl>>[c:1]O.Cl"),
("[c]Br", "[c:1]Br>>[c:1]O.Br"),
```

这些规则产生了化学上无效的单原子替换路线。当前 `plan_synthesis_recursive` 的 trivial 检测只检查 `reactants[0] == smiles`，但原子替换后的 reactant SMILES 与原始 SMILES 不同，导致漏检。

**影响维度**: 可合成性（route quality）
- 2/10 分子有 trivial 路线
- 复合评分中 `0.15 * route_quality` 给 trivial 分子 +0.10 优势分
- 挤占了真正有合成价值分子的 top-10 位置

---

## 4. 假设提出

### H031 — 单原子替换 Trivial 路线检测

```
假设ID: H031
瓶颈: 卤素交换规则产生的 OH↔Cl/Br 单原子替换路线未被识别为 trivial
文献支撑: LARC (Baker et al., 2025) — Agent-as-a-Judge 路线质量评审;
         标准药物化学实践 — 单官能团转化不是有效逆合成路线

───────────────── 推理链 ─────────────────
步骤 | 内容                                    | 置信度 | 推理方式 | 依据来源
S1   | 卤素交换规则仅做单原子替换，非真正合成路线 | 高     | 演绎     | 化学反应常识
S2   | 单原子替换的特征: 反应物与产物碳骨架相同  | 高     | 演绎     | RDKit 元素分析
     | 仅一个杂原子类型不同                     |        |          |
S3   | 检测此特征可标记为 trivial               | 高     | 演绎     | 源码分析
S4   | Trivial 标记 → route_quality=0           | 高     | 演绎     | score_route_quality 逻辑
     | → 复合评分中 route 权重归零              |        |          |
综合置信度 = 高

───────────────── 验证标准 ─────────────────
Q1 如果核心指标提升 < 5%，是否仍保留？
答：是
理由：修复 trivial 检测是正确性修复，不依赖于 BE 提升

Q2 如果指标下降，最可能的原因是什么？
答：Trivial 分子被移除后，top-10 被填补的分子可能 BE 偏低

Q3 本假设的最低可接受结果是什么？
答：Trivial ratio 从当前的实际 2/10 降至 0/10（知识库中声称的状态）

───────────────── 改进方案 ─────────────────
改动文件: src/synthesis_v2.py
改动内容:
  1. 新增函数 _is_single_atom_swap(smiles, reactants): 
     检测主反应物与原始分子是否仅差一个原子类型
  2. 在 plan_synthesis_recursive 中:
     当检测到单原子替换时，标记为 trivial
验证指标: run_pipeline 50 分子，对比 trivial ratio
```
