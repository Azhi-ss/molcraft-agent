# Competition Notes

## My Latest Submission Score (2026-05-19 17:34:49)
- score: 0.328537
- mol_score: 0.306887
- route_score: 0.379054
- sa_score: 0.65607
- validity_score: 1.0
- binding_score: 0.1766
- route_validity_score: 0.5
- starting_material_availability_score: 0.9
- sample_count: 10
- llm_score: 1.0

### Score Analysis
- mol_score = 0.8*0.1766 + 0.1*1.0 + 0.1*0.65607 = 0.14128 + 0.1 + 0.065607 = 0.306887 ✓
- route_validity_score=0.5 means 5/10 routes passed, 5/10 routes were zeroed
- sa_score=0.65607 means most molecules have SA>4 (only 2/10 have SA<4)
- binding_score=0.1766 is the biggest gap vs first place (0.8287)
- Validity is OK (1.0), starting materials OK (0.9), route_validity is the main route issue

## Gap Analysis vs First Place
| Dimension | Mine | First Place | Gap | Priority |
|-----------|------|-------------|-----|----------|
| binding_score | 0.1766 | 0.8287 | -0.6521 | 🔴 Critical (0.8 weight in mol_score) |
| sa_score | 0.6561 | 0.3034 | +0.3527 | 🟡 Mine higher but irrelevant (most mol SA>4 anyway) |
| validity_score | 1.0 | 1.0 | 0 | ✅ Both pass |
| route_validity_score | 0.5 | 1.0 | -0.5 | 🔴 High (half routes zeroed) |
| starting_material_availability | 0.9 | 1.0 | -0.1 | 🟡 Minor |
| llm_score | 1.0 | 0.5 | +0.5 | ✅ Mine better |
| mol_score | 0.3069 | 0.7933 | -0.4864 | Driven by binding_score gap |
| route_score | 0.3791 | 0.9497 | -0.5706 | Driven by route_validity gap |
| total score | 0.3285 | 0.8402 | -0.5117 | Two main gaps: binding + route_validity |

### Root Cause
1. binding_score gap (0.1766 vs 0.8287): Our Vina scores (~-9.9 best) likely get compressed by the normalization. First place likely has much higher raw Vina scores or better normalization alignment.
2. route_validity_score gap (0.5 vs 1.0): 5 of our 10 routes are invalid (final product mismatch / element imbalance / trivial A→A), causing route_score to be halved.

### Route Validity Detailed Diagnosis
5 routes zeroed due to balance_score = 0 (元素不守恒):
- #1 step2: Suzuki反应丢失 Br+B+2O (4 heavy atoms)
- #2 step2: Suzuki反应丢失 Br+B+2O (4 heavy atoms)
- #3 step1: Suzuki反应丢失 Br+B+2O (4 heavy atoms)
- #6 step1: Suzuki反应丢失 Br+B+2O (4 heavy atoms)
- #7 step3: Suzuki反应丢失 Br+B+2O (4 heavy atoms)
- #9 step1: Suzuki反应丢失 Br+B+2O (4 heavy atoms)

1 route zeroed due to trivial (A→A):
- #4: TRIVIAL route (product = reactant)

Pattern: All 6 invalid Suzuki steps write boronic acid as `OB(O)c1ccccc1` but
the balance check sees the B and O atoms from the boronic acid as "lost" because
Suzuki eliminates B(OH)2 + Br as waste, not tracked in the product. The scoring
system checks atom conservation strictly — waste/byproduct atoms must still be
accountable in the product side of the reaction SMILES.
- First place (from provided leaderboard screenshot):
  - score: 0.8402
  - mol_score: 0.7933
  - route_score: 0.9497
  - sa_score: 0.3034
  - validity_score: 1.0000
  - binding_score: 0.8287
  - route_validity_score: 1.0000
  - starting_material_availability_score: 1.0000
  - sample_count: 1.00
  - llm_score: 0.50
  - date/time: not visible in screenshot

## Scoring Rules (FAQ Summary)
- Molecular score (weight 0.7):
  - 0.8 × binding_score (AutoDock Vina)
  - 0.1 × validity_score (0 or 1)
  - 0.1 × sa_score (SAScore > 4 → 0; < 4 lower is better)
- Route score (weight 0.3):
  - 0.55 × route_validity_score
  - 0.30 × starting_material_availability_score (DB hit; else SAScore fallback)
  - 0.05 × step_penalty_score
  - 0.05 × convergence_score
  - 0.05 × balance_score

## Zeroing Conditions (FAQ)
- If validity_score = 0, molecular score becomes 0.
- If route final product ≠ designed molecule, route score becomes 0.
- If balance_score = 0, route score becomes 0.

## FAQ: route_validity_score = 1.0 but route_score = 0
Possible causes:
- Final product is not the designed molecule.
- Any step introduces elements not present in reactants (element appears out of nowhere).
- Any step has product identical to reactant (A → A).

---

## Bug Fix: Suzuki Reaction Atom Imbalance (H019)

**Status**: ✅ FIXED  
**Before**: route_validity_score = 0.5 (5/10 routes failed)  
**Expected After Fix**: route_validity_score ≈ 1.0

### Problem
Suzuki coupling逆合成反应的反应SMILES不包含副产物，导致原子守恒检查失败：
```
Ar-Br + Ar'-B(OH)2 → Ar-Ar'
反应物原子数 = Ar + Br + Ar' + B + 2O  
产物原子数 = Ar + Ar'
差值 = 4 个原子 ❌
```

### Root Cause
Suzuki反应实际生成副产物 `B(OH)3 + HBr`，但逆合成SMILES没有在产物侧写出这些副产物。评分系统的`balance_score`严格检查原子守恒，差值>2即归零。

### Fix Applied (`src/synthesis_v2.py`)
1. 在`_try_route`函数中自动检测反应类型并添加对应的副产物
2. Suzuki反应：产物侧添加`OB(O)O + Br`平衡4个重原子
3. Friedländer喹啉合成：产物侧添加`H2O + H2O`平衡脱水反应

### Verification
- 修复前: 1/9规则不平衡
- 修复后: 0/9规则不平衡 ✓
- 联芳基类分子测试全部通过原子平衡检查

### Leaderboard Analysis (2026-05-19)

### Key Finding: sample_count = 1.00 for ALL teams

| Team | score | mol_score | route_score | binding_score | sa_score | sample_count |
|------|-------|-----------|-------------|---------------|----------|--------------|
| 1st | 0.8402 | 0.7933 | 0.9497 | 0.8287 | 0.3034 | 1.00 |
| 2nd | 0.7896 | 0.6997 | 0.9895 | 0.7196 | 0.3034 | 1.00 |
| 3rd | 0.7210 | 0.5864 | 0.9900 | 0.5775 | 0.3034 | 1.00 |
| 4th | 0.7000 | 0.5632 | 0.9806 | 0.5538 | 0.3034 | 1.00 |
| 5th | 0.6272 | 0.4699 | 0.9716 | 0.4248 | 0.3034 | 1.00 |
| 6th | 0.5941 | 0.4218 | 0.9861 | 0.3597 | 0.3034 | 1.00 |
| 7th | 0.5681 | 0.3853 | 0.9877 | 0.3067 | 0.3034 | 1.00 |
| 8th | 0.5600 | 0.3805 | 0.9805 | 0.3006 | 0.3034 | 1.00 |
| **Ours** | **0.3285** | **0.3069** | **0.3791** | **0.1766** | **0.6561** | **10** |

### Verified Formulas
- `mol_score = 0.8*binding + 0.1*validity + 0.1*sa` ✓ matches all 8 teams (diff < 0.0001)
- `total_score = 0.7*mol + 0.3*route` ✓ matches all 8 teams

### Critical Observations
1. **sample_count=1.00 for ALL teams**: Every team on the leaderboard has exactly 1 molecule evaluated. Our submission has sample_count=10.
2. **All teams have route_validity_score=1.0, starting_material_availability=1.0**: Route quality is perfect for all ranked teams.
3. **All teams have validity_score=1.0, llm_score=0.5**: Validity is baseline; LLM score varies (we have 1.0, they have 0.5).
4. **sa_score is bimodal**: Most teams have 0.3034 (SA > 4 → low score), ours is 0.6561 (SA < 4 → higher score). But since SA only has 0.1 weight in mol_score, this difference is minor.
5. **binding_score is the ONLY differentiator**: Teams are ranked almost exactly by binding_score. The gap from 8th place (0.3006) to 1st place (0.8287) is huge.

### Hypothesis: Competition Evaluates 1 Molecule Only

If the competition only scores 1 molecule per submission:
- Our `binding_score=0.1766` might correspond to our best molecule's Vina in the competition's scoring system
- Our local best Vina = -9.335. If competition uses `range=15`: score = 9.335/15 = 0.622
- But our actual score is 0.1766, suggesting either:
  a) Our local Vina is systematically different from competition Vina (different docking params)
  b) The system doesn't pick our best molecule
  c) range is much larger than 15 (our calibration suggests ~48.6)

### Strategic Implications
1. **Submit 1 molecule vs 10**: If sample_count=1 is real, submitting only our best molecule could change scoring. Need to test.
2. **binding_score gap is THE priority**: 0.1766 vs 0.8287 = -0.6521 gap. At 0.8 weight in mol_score and 0.7 weight in total, this alone accounts for -0.365 total score gap.
3. **Route validity is solved**: Suzuki fix should bring route_validity from 0.5 to 1.0.
4. **SA score is a red herring**: Our SA is "better" (0.6561 vs 0.3034) but this is actually BAD — it means our molecules are TOO SIMPLE (low SA = easy to synthesize = maybe less drug-like). But SA only has 0.1 weight.

---

## Next Steps
1. ✅ Suzuki fix verified working — regenerate result.csv with fixed routes
2. 🔄 Test sample_count=1 hypothesis: submit single best molecule
3. 🔄 Continue investigating binding_score gap (local Vina vs competition Vina)
4. 🔄 Generate higher-binding molecules (current best: -9.335, need ~-12+ to compete)

---

## Local Scoring System — Calibration Results

### Baseline (default params: range=15, SA step cutoff=4)
| Metric | Local | Actual | Delta |
|--------|-------|--------|-------|
| binding_score | 0.572 | 0.1766 | +0.395 (3.2x over!) |
| sa_score | 0.043 | 0.656 | -0.613 |
| total_score | 0.667 | 0.329 | +0.338 |

### Calibrated (range=48.59, SA inverted scale=10)
| Metric | Local | Actual | Delta |
|--------|-------|--------|-------|
| binding_score | 0.1766 | 0.1766 | 0 ✅ |
| sa_score | 0.564 | 0.656 | -0.092 🟡 |
| total_score | 0.482 | 0.329 | +0.153 |

### Key Findings
1. **binding_score 归一化**：`clipped_linear` with `threshold=0, range≈48.6`。range=15 严重低估了归一化分母，比赛系统的 range 大约是 48-50。这意味着比赛把 Vina -8.5 映射到仅 0.17 分，非常压缩。
2. **sa_score 归一化**：不是简单的 step 函数（SA>4→0）。比赛用的是某种 inverted 线性函数，但我们的本地 SA 估计器（`estimate_sa_score`）与比赛的 SAscore 可能不同。需要进一步校准。
3. **route_validity_score**：本地预测 0.9（Suzuki 修复后），实际 0.5（修复前）。下次提交应验证修复效果。
