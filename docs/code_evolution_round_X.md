# 代码演进 — Round X (新会话 H031)

## 日期: 2026-05-21

---

## 改动概述

**假设ID**: H031 — 单原子替换 Trivial 路线检测

**改动文件**: `src/synthesis_v2.py`

---

## 具体修改

### 1. 新增 `_is_single_atom_swap()` 函数 (第 633-695 行)

检测逆合成反应是否仅为化学上 trivial 的单原子替换（如 OH↔Cl/Br 交换）。

**算法**:
- 从反应物列表中找出重原子数最多的「主反应物」（排除 Cl、Br、F 等小试剂）
- 比较原始分子与主反应物的重原子数（必须相同）
- 统计元素组成差异（排除 H）
- 若恰好两种元素各相差 ±1 且净差为 0 → 单原子替换

**示例**:
- `ClC6H4-C(=O)-...` → 主反应物 `HOC6H4-C(=O)-...` + `Cl` → 检测为 swap (Cl↔O)
- `c1ccccc1-c2ccccc2` → `c1ccccc1Br` + `B(OH)O-c2ccccc2` → 非 swap（碳原子数不同）

### 2. 修改 `plan_synthesis_recursive` trivial 判定 (第 757-759 行)

**修改前**:
```python
is_trivial = all_trivial and len(reactants) == 1 and reactants[0] == smiles
```

**修改后**:
```python
is_trivial = all_trivial and len(reactants) == 1 and reactants[0] == smiles
if not is_trivial and all_trivial:
    is_trivial = _is_single_atom_swap(smiles, reactants)
```

**逻辑解释**:
- 保留原有 exact-match trivial 检测
- 新增: 当所有子路线都是 trivial（all_trivial=True）且检测到单原子替换 → 整体标记为 trivial
- 不影响正常多步路线的检测

---

## 文献支撑

- LARC (Baker et al., 2025): Agent-as-a-Judge 路线质量评审
- 标准药物化学实践: 单官能团转化（OH→Cl）不是有效逆合成路线

---

## 编译验证

```bash
$ python3 -m py_compile src/synthesis_v2.py
# 通过，无语法错误
```

## 功能验证

| 测试用例 | 修改前 | 修改后 | 期望 |
|---------|--------|--------|------|
| `O=C1c2ccccc2NNc2cccc(Cl)c21` | trivial=False | trivial=True ✅ | True |
| `CC1CCC2(CC1)Cc1cccc(Cl)c1C2` | trivial=False | trivial=True ✅ | True |
| `Nc1ccccc1-c1cc...cc1O` (Suzuki) | trivial=False | trivial=False ✅ | False |
