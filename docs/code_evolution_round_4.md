# 代码演进报告 Round 4 — H004

**日期:** 2026-05-18
**假设:** H004 适应度比例选择 + 精英保留
**状态:** ❌ REJECTED — 回滚建议，但先保留供后续参考

---

## 修改文件

### `tools/pipeline.py`

**变更 1: `_generate_offspring()` — 适应度比例选择 (Lines 310-358)**

```diff
- seed = random.choice(seeds)
+ # H004: 计算适应度权重（基于结合能）
+ energies = [s.get("binding_energy", 0.0) for s in seeds]
+ fitness = [-e for e in energies]
+ min_fit = min(fitness) if fitness else 0.0
+ shifted = [f - min_fit + 0.1 for f in fitness]
+ total = sum(shifted)
+ weights = [s / total for s in shifted] if total > 0 else None
+ seed = random.choices(seeds, weights=weights, k=1)[0]
```

**变更 2: 主循环 — 精英保留 (Lines 182-186)**

```diff
- n_seeds = min(20, len(successful))
- current_seeds = successful[:n_seeds]
+ N_ELITE = min(5, max(1, len(successful) // 10))
+ elite = successful[:N_ELITE]
+ current_seeds = list(elite) + successful[N_ELITE:n_seeds]
```

---

## 回滚方案

如需回滚，将 `tools/pipeline.py` 中上述两处恢复为原始代码。变更已用 `# H004:` 标注，易于定位。
