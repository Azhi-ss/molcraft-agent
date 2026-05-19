# 代码演进报告 — Round 8 (H016)

## 假设ID
H016: 添加 Diels-Alder 环加成逆反应规则

## 修改文件
`src/synthesis_v2.py`

## 修改内容摘要

在 `RETRO_RULES` 列表中添加了 3 条 Diels-Alder 环加成逆反应规则：

### 1. 环己烯 Diels-Alder 逆反应（通用型）
```python
("[C;R;!a]1=[C;R;!a][C;R;!a][C;R;!a][C;R;!a][C;R;!a]1",
 "[C:1]1=[C:2][C:3][C:4][C:5][C:6]1>>[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]"),
```
- 匹配: 六元非芳香碳环单烯烃
- 产物: 1,3-丁二烯片段 + 乙烯片段
- 验证: tetrahydrophthalimide → butadiene + maleimide ✓
- 验证: 4-phenylcyclohexene → 1-phenylbutadiene + ethene ✓

### 2. 二氢吡喃 Hetero-Diels-Alder 逆反应（变体A）
```python
("[O,S;R]1[C;R;!a]=[C;R;!a][C;R;!a][C;R;!a][C;R;!a]1",
 "[O,S:1]1[C:2]=[C:3][C:4][C:5][C:6]1>>[O,S:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]"),
```
- 匹配: 双键紧邻杂原子的二氢吡喃/二氢噻喃
- 产物: α,β-不饱和羰基 + 乙烯
- 验证: 3,4-dihydro-2H-pyran → acrolein + ethene ✓

### 3. 二氢吡喃 Hetero-Diels-Alder 逆反应（变体B）
```python
("[O,S;R]1[C;R;!a][C;R;!a][C;R;!a]=[C;R;!a][C;R;!a]1",
 "[O,S:1]1[C:2][C:3][C:4]=[C:5][C:6]1>>[O,S:1][C:2]=[C:3][C:4]=[C:5].[C:6]"),
```
- 匹配: 双键远离杂原子的二氢吡喃/二氢噻喃
- 覆盖双键在不同位置的环系变体

## 文献依据
- **Deep Lead Optimization (JACS 2024)**: 稠环体系需要专门断键策略，Diels-Alder 是合成六元碳环/桥环的核心反应
- **LARC (Baker et al., 2025)**: 规则覆盖率决定逆合成质量
- **综述第3.2节**: 周环反应是药物化学关键转化

## 设计注意事项
- 移除了 `!$(C=*)` 约束（过于严格，阻止匹配烯烃碳）
- 杂原子规则分两个变体覆盖双键不同位置
- 跳过了二氢吡啶酮和稠合环己烯规则（SMARTS 复杂度高，留待后续迭代）
- 原子守恒 (6→4+2) 通过 H014 mass balance 自动验证

## RETRO_RULES 数量变化
- 修改前: 52 条规则
- 修改后: 55 条规则 (+3)

## 已知限制
- `_is_simple_molecule` 在 ≤10 原子时提前返回，阻止对小分子的 DA 断键
- 稠合环己烯（decalin 类）和桥环体系（norbornene）的精确 SMARTS 仍需后续迭代
- 杂原子规则变体B的产品 SMILES 可能需要进一步验证
