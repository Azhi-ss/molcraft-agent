---
name: molcraft-synthesis-rules
description: 逆合成 SMARTS 规则设计指南。当需要扩充 synthesis_v2.py 的 REACTION_RULES 时读取。包含 RDKit ReactionFromSmarts 语法、断键启发式、新规则设计步骤。
---

# 逆合成 SMARTS 规则设计指南

扩充 synthesis_v2.py 的 `REACTION_RULES` 列表时使用。

## RDKit 语法

```python
from rdkit.Chem import AllChem
# 格式: '反应物>>产物'
# [c:1] 标记原子映射，序号对应
rxn = AllChem.ReactionFromSmarts('[c:1][c:2]>>[c:1]Br.[c:2]B(O)O')  # Suzuki 逆反应
# REACTION_RULES 中每条: (匹配子结构SMARTS, 反应SMARTS)
```

## 断键启发式

| 键类型 | 可用的逆合成反应 |
|--------|----------------|
| 新 C-C | aldol, Claisen, Michael, Wittig, Grignard, Diels-Alder, Suzuki |
| 新 C-O | 酯化, 醚化(Williamson), Mitsunobu |
| 新 C-N | 还原胺化, 酰胺偶联, Buchwald, SNAr, Gabriel |
| 稠环断开 | Diels-Alder 逆反应, Friedel-Crafts 环化逆反应 |
| 桥环断开 | 分子内 SN2 逆反应, 自由基环化逆反应 |
| 螺环断开 | pinacol 重排逆反应, 半缩酮逆反应 |
| 杂环(吡啶/嘧啶/吲哚) | 逆 Pictet-Spengler, 逆 Bischler-Napieralski |

## 新规则设计步骤

1. 分析当前 result.csv/result.log 中 trivial route 分子的结构
2. 看哪些键类型被 BRICS 回退了
3. 对照上面启发式判断正确断键方式
4. 用 RDKit SMARTS 写匹配子结构和逆合成反应
5. 加到 synthesis_v2.py 的 `REACTION_RULES`: `("子结构SMARTS", "反应SMARTS")`
6. 跑 pipeline 验证 trivial 比例变化

## 关键陷阱

- SMARTS 中 `[c;R][c;R]` 会匹配环内 C-C 键（如苯环上的相邻碳）→ 用 `!@` 确保只匹配环间键: `[c;R]!@[c;R]`
- 产物 SMILES 中的原子映射编号必须连续（从 1 开始）
- 新规则可能意外匹配已有分子的其他部分 → 在子结构 SMARTS 中加环境约束（如 `[NH1]` 而非 `[N]`）
- 加规则后检查 _is_simple_molecule 是否误判复杂分子为简单