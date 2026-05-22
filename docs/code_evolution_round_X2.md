# 代码演进 — Round 2 (H032)

## 日期: 2026-05-21

---

## 改动概述

**假设ID**: H032 — Post-Synthesis Validity Filter

**改动文件**: `tools/pipeline.py`

---

## 具体修改

### 后合成有效性过滤器 (第 296-312 行)

**修改位置**: 在复合评分排序后、top-N 选择前插入过滤器。

**新增代码**:
```python
# H032: Post-synthesis validity filter
# Exclude candidates with smiles>>smiles trivial routes
valid_candidates = [c for c in scored_candidates 
                    if c["route"] != f"{c['mol_smiles']}>>{c['mol_smiles']}"]
n_filtered = len(scored_candidates) - len(valid_candidates)
if n_filtered > 0:
    log(f"H032 后过滤: 排除 {n_filtered} 个无有效合成路线的分子")
if len(valid_candidates) >= n_top:
    final_top = valid_candidates[:n_top]
else:
    final_top = valid_candidates
```

**逻辑**:
1. 从已排序的候选池中排除 `route == smiles>>smiles` 的分子
2. 有效候选充足时取 top-N
3. 有效候选不足时保留全部（不降级填充无效分子）

---

## 文献支撑

- LARC (Baker et al., 2025): Agent-as-a-Judge — 无有效逆合成路线的分子应被排除
- MOOSE-Chem (Yang et al., 2025): 多目标优化需满足所有约束条件

---

## 编译验证

```bash
$ python3 -m py_compile tools/pipeline.py
# 通过 ✅
```
