# 诊断报告 — Round X (H029 → H030)

## 1. 当前状态

- **H029 VERIFIED**: 多构象对接增强，Best BE -10.099 (+9.9%), Avg BE -9.412 (+11.8%)
- **副作用**: Trivial ratio 0/10 → 2/10 ⚠️
- **根因**: Vina 疏水偏置导致多环芳烃泛滥 → 合成规则无法覆盖 → 逆合成退化为 trivial

## 2. 外部知识搜索

### 搜索关键词
- "Vina docking hydrophobic bias scoring correction"
- "LogP penalized composite scoring drug discovery"
- "AutoDock Vina scoring function limitations kinase"

### 发现
1. **Vina Scoring Biases**（molcraft-vina-strategies skill）:
   - Hydrophobic overestimation: 过度奖励非极性表面积埋藏
   - Desolvation penalty weak: 隐式溶剂模型低估极性基团去溶剂化代价
   - Entropic penalty flat: 配体熵惩罚约为每个可旋转键 ~0.3 kcal/mol

2. **Open-Source Molecular Docking Review** (IJMS, 2026):
   - Vina 家族对接流程需要"注释氢键、疏水接触、盐桥"
   - 化学偏置的基准集可能显著影响结果

3. **Multi-parameter optimization of kinase inhibitors** (Kabiri, 2025):
   - 利用物理化学方法进行多参数优化
   - 铰链区结合和疏水后袋是评分函数优化的关键维度

## 3. 源码分析

### 关键排序路径

```
generate_with_docking_guidance:
  L1141: docked_batch.sort(key=BE)          ← 纯BE种子选择
  L1207: unique.sort(key=BE)               ← 最终池纯BE排序

run_evolutionary_pipeline:
  L174:  unique_docked.sort(key=BE)         ← 纯BE排序
  L199:  consensus_results.sort(key=BE)     ← 纯BE共识排序
  L261:  composite = 0.8*BE_norm + 0.2*route_quality  ← H012复合（但trivial route_quality=0）
```

### 问题分析

分子A（trivial, BE=-10.0, route_quality=0.0）:
- BE_norm ≈ 1.0, composite = 0.8×1.0 + 0.2×0.0 = **0.80**

分子B（non-trivial, BE=-9.2, route_quality=0.8）:
- BE_norm ≈ 0.75, composite = 0.8×0.75 + 0.2×0.8 = **0.76**

→ 即使有完美合成路线，BE 差距 0.8 kcal/mol 就足以让 trivial 分子胜出！

### LogP 分布推测

Vina 偏好的多环芳烃 LogP 通常在 3.5-5.0 范围，而可合成的铰链结合分子 LogP 通常在 1.5-3.5。LogP 惩罚可以直接针对这一差距。

## 4. 假设提出

```
假设ID: H030
瓶颈: Vina 疏水偏置导致 LogP>4 的多环芳烃主导 top-N，其逆合成因规则覆盖不足退化为 trivial route（当前 2/10）
文献支撑: Vina Scoring Known Biases (molcraft-vina-strategies skill); Lipinski Rule of 5; MOOSE-Chem multi-objective optimization

───────────────── 推理链 ─────────────────
步骤 | 内容                                                    | 置信度 | 推理方式 | 依据来源
S1   | Vina 过度奖励疏水表面积埋藏，LogP>4 分子被高估 0.5-1.0 kcal/mol | 高   | 文献    | Vina skill §Scoring Function Known Biases
S2   | 当前复合评分 0.8×BE+0.2×route 无法抵消 BE 偏置            | 高   | 演绎    | 源码分析（见上）
S3   | 添加 LogP 惩罚项可降低过度疏水分子的复合评分              | 高   | 文献    | Lipinski 规则; 多参数优化 (Kabiri, 2025)
S4   | 将 LogP penalty 权重设为 10%，BE 降为 75% 不影响主导      | 中   | 演绎    | 经验权重设计，需实验验证
综合置信度 = 中（权重比例需实验调优，S4 为经验推断）

───────────────── 验证标准 ─────────────────
Q1 如果核心指标提升 < 5%，是否仍保留？
答：是
理由：即使 BE 不变，trivial ratio 从 2/10 降至 0-1/10 就是净收益。BE 可能因疏水分子被降权而轻微下降，但预期可合成性改善。

Q2 如果指标下降，最可能的原因是什么？
答：LogP penalty 权重过高，排除了合理疏水的铰链结合分子（LogP 3-4 的 kinase hinge binder 是正常的）
理由：过度惩罚可能排除药效团必需的芳香体系

Q3 本假设的最低可接受结果是什么？
答：Trivial ratio ≤ 1/10，Best BE ≥ -9.5 kcal/mol（即 trivial 改善同时 BE 不低于基线 94%）

───────────────── 改进方案 ─────────────────
改动文件: tools/pipeline.py
改动内容: 
  1. 添加 compute_logp_score(logp) 函数（LogP≤2→1.0, LogP≥5→0.0, 线性衰减）
  2. 修改 H012 复合评分公式: 0.75*BE_norm + 0.15*route_quality + 0.10*logp_score
  3. 确保 scored_candidates 保留 logp 字段
验证指标: pipeline 跑 n_generate=50, n_top=10, strategy=mutate, n_generations=2
         对比 H029 基线：Best BE -10.099, Avg -9.412, Trivial 2/10
```
