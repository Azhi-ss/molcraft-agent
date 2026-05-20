# Bottleneck Diagnosis — Round 1 (TYK2 5C01 Corrected Coordinates)

## Baseline Results (Corrected Docking)

| Metric | Value |
|--------|-------|
| Best BE | -10.119 kcal/mol |
| Avg BE | -9.417 kcal/mol |
| Trivial ratio | 0/10 |
| Strategy | mutate, docking_guidance=on |

## Molecular Analysis

Top 10 molecules:
1. Diaryl benzophenone (BE=-10.119) — Suzuki + Friedel-Crafts, 2-step
2. Imino-biaryl benzophenone (BE=-9.762) — Imine + Suzuki
3. Amino-biaryl benzophenone (BE=-9.725) — Suzuki + F-C
4. Styrenyl biaryl (BE=-9.704) — Wittig + Suzuki
5. Biaryl benzophenone (BE=-9.537) — F-C + Suzuki
6. Amino-biaryl benzophenone (BE=-9.509) — F-C + Suzuki
7. Aminopyridine-biaryl (BE=-9.351) — Suzuki + F-C
8. Sulfonamide-biaryl (BE=-9.085) — Amide + Suzuki
9. Thienyl-biaryl (BE=-9.027) — F-C + Suzuki
10. Benzoxazole-phenyl (BE=-8.354) — Halogenation + Suzuki

## Bottleneck Identified: Scaffold Homogeneity

**Observation**: 9/10 molecules are benzophenone/biaryl chemotypes. Despite the SCAFFOLDS library containing 90+ entries including kinase-specific heterocycles (purines, pyrazolopyrimidines, azaindoles, imidazopyridines), none appeared in the top results.

**Root Cause**: `random.choice(seeds)` in `generate_molecules()` treats all scaffolds uniformly. With 90+ scaffolds, kinase-privileged entries (~8 entries: purine, azaindole, pyrazolopyrimidine, imidazopyridine, pyrazolopyrimidine, indolizine, quinazoline, benzimidazole) each have only ~1% selection probability per seed. The benzophenone/biaryl scaffolds produce better initial docking scores, so they dominate the evolutionary process.

**Impact**: 
- Limited chemical diversity in output
- Misses TYK2 hinge-binding pharmacophore (NH...O=C hydrogen bonds with hinge residues Glu979/Met981)
- All molecules bind through similar interactions (mostly hydrophobic + polar contacts in the DFG-out pocket)

---

## Proposed Hypothesis: H021

```
假设ID: H021
瓶颈: 分子骨架同质化——90%+输出为二苯甲酮/联芳基类型，激酶铰链结合杂环骨架完全缺失
文献支撑: 
  - Deep Lead Optimization (JACS 2024, §3.1 Scaffold Hopping): 
    "替换核心骨架同时保留有利取代基"是先导化合物优化的核心子任务
  - MOOSE-Chem (Yang et al., 2025): 
    "Diverse initial population is essential for evolutionary search to avoid premature convergence"
  - Coscientist (Boiko et al., 2023):
    正确的化学空间导航需要靶点知识指导骨架选择

───────────────── 推理链 ─────────────────
步骤 | 内容                                  | 置信度 | 推理方式 | 依据来源
S1   | TYK2 是激酶，铰链结合是经典药效团   | 高     | 文献     | TYK2 晶体结构文献
S2   | 激酶铰链结合杂环（嘌呤/氮杂吲哚等）  | 高     | 文献     | 激酶抑制剂化学综述
      | 在 SCAFFOLDS 中已存在但概率稀释     |         |          |
S3   | 对激酶骨架加权采样会增加其出现概率   | 高     | 演绎     | 从 S1,S2
S4   | 增加激酶骨架多样性可能发现更高        | 中     | 演绎     | 从 S2,S3
      | 亲和力的铰链结合分子                   |         |          |
S5   | 加权采样不破坏现有成功路径            | 中     | 演绎     | 从 S3

综合置信度 = 中
S4 为"中"——激酶骨架是否在 TYK2 口袋中真的有更好亲和力，需实验验证
S5 为"中"——过度加权可能挤出已验证的高分路径

───────────────── 验证标准 ─────────────────
Q1 如果核心指标提升 < 5%，是否仍保留？
答：是
理由：即使 BE 不提升，化学多样性改善本身有价值（探索新的化学空间）

Q2 如果指标下降，最可能的原因是什么？
答：激酶杂环骨架在本口袋中亲和力不如二苯甲酮类型
理由：TYK2 的 ATP 口袋可能对此类骨架不敏感，或铰链区构象不适合经典铰链结合

Q3 本假设的最低可接受结果是什么？
答：top10 中出现 ≥2 个含激酶杂环骨架的分子，且 avg BE > -8.5
（允许 BE 轻微下降换取化学多样性）

───────────────── 改进方案 ─────────────────
改动文件: src/generator.py
改动内容: 
  1. 将 SCAFFOLDS 分为 KINASE_SCAFFOLDS (激酶铰链结合杂环) 和 GENERAL_SCAFFOLDS 两层
  2. 在 generate_molecules() 种子选择中，50% 概率从 KINASE 层采样
  3. 维持所有其他变异算子不变
验证指标: pipeline 跑 50 分子，对比改前改后的激酶骨架出现率 + 结合能 + trivial 比例
```

---

## Decision

Move to Stage 3 (Code Evolution) to implement H021.
