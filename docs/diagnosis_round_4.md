# 瓶颈诊断报告 Round 4

**诊断日期:** 2026-05-18
**基于:** literature_analysis_round_4.md + 源码审查 + 历史实验数据

---

## 瓶颈 H004: 进化选择缺乏适应度压力

### 症状

| 指标 | 数值 | 诊断 |
|------|------|------|
| 化学多样性 | 80% 分子为磺酰胺-双芳环模式 | 选择压力不足，搜索未充分探索空间 |
| 全局最佳结合能 | -8.56 kcal/mol | 三迭代无明显提升趋势（缺乏方向性进化） |
| Top-10 平均 | -8.06 kcal/mol | 退化方差小，说明局部极值饱和 |

### 追踪到代码

**文件:** `tools/pipeline.py`  
**函数:** `_generate_offspring(seeds, n_offspring_per_seed)`  

**问题代码 (Line 308-337):**

```python
def _generate_offspring(seeds, n_offspring_per_seed):
    # ...
    while len(offspring) < len(seeds) * n_offspring_per_seed:
        seed = random.choice(seeds)          # ← LINE 323: 均匀随机！
        # 所有种子被选中的概率相等，无论结合能高低
```

**关键缺陷:** `random.choice(seeds)` 对所有种子施加**零选择压力**。结合能 -8.56 kcal/mol 的种子与 -7.0 的种子有完全相同的概率被选中变异。这在进化算法中属于根本性设计错误——**选择压力为零，等于无选择的随机行走**。

### 为什么造成了瓶颈

1. **无梯度利用:** top 种子（优秀分子）没有比 bottom 种子（平庸分子）获得更多复制机会
2. **搜索效率浪费:** 约 50% 的变异尝试来自结合能较差的种子，等于在向错误方向探索
3. **收敛假象:** 8/10 分子收敛于磺酰胺模式，不是因为磺酰胺最优，而是因为早期随机漂移导致种群同质化，之后的选择压力不足以突破
4. **缺失精英保留:** 每代最佳分子不会被直接保留到下一代（line 159: mols 完全由变异产生），存在退化风险

### 文献支撑

| 文献 | 相关陈述 |
|------|---------|
| MOOSE-Chem (Yang et al., 2025) | "fitness function evaluates the quality of the hypothesis to guide selection toward promising regions" |
| Survey §3.2 | "Diverse initial population is essential for evolutionary search" + "fitness-guided selection" |
| MolLEO (Wang et al., 2024b) | "selection pressure is a key driver of evolutionary efficiency" |

### 拟议修复

**范围:** `tools/pipeline.py` 单文件，两个修改点:

**修改1:** `_generate_offspring()` — 适应度比例选择 (line 323)

当前:
```python
seed = random.choice(seeds)
```

替换为:
```python
# 适应度比例选择: 结合能越低（越负），被选中概率越高
# 使用 rank-based 或 softmax over binding_energy
energies = [s.get("binding_energy", 0) for s in seeds]
# 转换为最大化问题: neg_energy = -binding_energy
neg_energies = [-e for e in energies]
# 平移使最小值为 0.1（避免零权重）
min_val = min(neg_energies)
shifted = [e - min_val + 0.1 for e in neg_energies]
weights = [e / sum(shifted) for e in shifted]
seed = random.choices(seeds, weights=weights, k=1)[0]
```

**修改2:** 主循环 (line 182-184) — 精英保留

当前:
```python
n_seeds = min(20, len(successful))
current_seeds = successful[:n_seeds]
```

替换为:
```python
n_seeds = min(20, len(successful))
# 精英保留: top N_ELITE 直接保留，其余从全部 successful 中按 fitness 选择
N_ELITE = min(5, max(1, len(successful) // 10))
elite = successful[:N_ELITE]
current_seeds = list(elite) + successful[N_ELITE:n_seeds]
```

### 预期效果

| 指标 | 当前基线 | 预期改善 | 理由 |
|------|---------|---------|------|
| 最佳结合能 | -8.56 | -8.7 ~ -9.0 | 优秀分子获得更多变异机会 |
| Top-10 平均 | -8.06 | -8.2 ~ -8.5 | 选择压力促进系统化改进 |
| 化学多样性 | 磺酰胺主导 | 改善 | 精英保留防止漂移 + 适应度选择引导探索 |
| Trivial route | 10% | ≤20% | 无直接关系，但更好的分子可能有更简单的合成 |

### 风险与局限

| 风险 | 缓解 |
|------|------|
| 过早收敛 | elite 仅保留 5 个，配合 rank-based 而非 hard-cutoff 选择保留探索性 |
| 权重映射对正值 binding_energy 不稳定 | 使用平移 + 最小值修正确保所有权重 > 0 |
