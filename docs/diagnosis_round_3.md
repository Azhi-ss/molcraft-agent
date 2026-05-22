# Diagnosis Round 3 — H040

## Stage: 诊断 (2026-05-22)

---

## 1. REJECTED Hypotheses Check

| ID | Reason Rejected | Avoid? |
|----|----------------|--------|
| H019 | Suzuki byproduct atom balance → trivial regression | ✅ Avoid |
| H025 | Polar substituents → Vina penalizes polarity | ✅ Avoid |
| H026 | Large aromatic scaffolds → too large for TYK2 pocket | ✅ Avoid |
| H036 | Larger batch/top_k → noise degraded BE | ✅ Avoid |

**Conclusion**: No conflict with proposed H040 direction.

---

## 2. External Knowledge Search

### LKM Claims Retrieved:

1. **TYK2 hinge Tyr vs JAK1 Phe** (gcn_e1c03592ed8942c1): TYK2/JAK2/JAK3 have Tyr at hinge position → polar OH can form ligand-specific H-bond. JAK1 has Phe → lipophilic hinge. This means TYK2 prefers polar hinge-binders.

2. **TYK2 JH2 pseudokinase domain druggability** (gcn_822926bd3ca94881): 5C01 JH2 domain binds pyrazine inhibitor at Kd=0.25μM — confirms ATP pocket is tractable for drug design.

3. **Sulfonamide as alternative hinge binder** (gcn_0136b5496bab4c01): Docking suggests sulfonamide moiety can occupy kinase hinge region, providing alternative to amide-based hinge binding.

4. **Multi-parameter optimization** (LARC, Baker et al., 2025): Synthesis tractability should be weighted alongside binding affinity in lead selection.

### arXiv: No TYK2-specific recent papers found.

---

## 3. Current State Analysis

| Metric | H037 (baseline) | H038 (n_conf=2) | Δ |
|--------|:---:|:---:|:---:|
| Best BE | -10.139 | -9.795 | -3.4% |
| Avg BE | -9.021 | -9.080 | +0.7% |
| Trivial | 2/10 | 0/10 | ✅ |

Chemistry in H038: Diverse — quinoline-pyrazolidine, sulfonamide-pyrrolopyridine, benzophenone-biaryl, quinoline-quinoline, tetrahydroquinoline-piperazine. Only 1-2/10 are sulfonamide-based. Good scaffold diversity.

---

## 4. Bottleneck Diagnosis

**Primary bottleneck**: The pipeline is approaching the Vina scoring ceiling for TYK2 5C01. Best BE oscillates between -9.5 and -10.1 across runs with different random seeds. Further raw BE improvements are diminishing returns.

**Secondary observation**: The current composite scoring (0.75×BE + 0.15×route + 0.10×logp) heavily weights Vina score (56% of total via 0.75×BE in composite × 0.7 mol in total). This may select molecules that Vina likes but have suboptimal synthesis routes.

**Opportunity**: Shift optimization target from "best BE" to "best overall drug candidate" — improving synthesis tractability and drug-likeness while maintaining competitive BE.

---

## 5. Hypothesis: H040 — Route-Quality-Weighted Composite Scoring

### Statement
Increase route quality weight in composite scoring from 0.15→0.25, decreasing BE weight from 0.75→0.65. This rewards molecules with better retrosynthetic routes (fewer steps, better atom economy, more accessible building blocks) at a modest cost to raw Vina score.

### New scoring: `0.65×BE_norm + 0.25×route_quality + 0.10×logp_score`

### Literature Basis
- LARC (Baker et al., 2025): Multi-parameter optimization including synthetic tractability is essential for practical drug discovery
- Coscientist (Boiko et al., 2023): "Success in drug discovery requires balancing multiple objectives"
- Route quality scoring (H012) is already implemented and validated — we just need to increase its weight

### Expected Effect
- Best BE: may decrease 0.2-0.5 kcal/mol (acceptable tradeoff)
- Avg route quality: should improve (fewer steps, better atom economy)
- Trivial ratio: should remain 0/10
- Overall: Better real-world drug candidates with more practical synthesis routes

### Risk Assessment
- Low risk: Only changes two coefficients in existing scoring formula
- If BE degrades > 1.0 kcal/mol → REJECT and revert
- If route quality improves and BE degradation < 0.5 → VERIFY

### Implementation
- Single line change in `tools/pipeline.py` line 311
- Change `0.75 * be_norm + 0.15 * c["route_quality"] + 0.10 * logp_score` → `0.65 * be_norm + 0.25 * c["route_quality"] + 0.10 * logp_score`

---

## 6. Tool Chain Review

- No external tool changes needed
- No new dependencies
- RDKit skill not needed (no new cheminformatics)
- Vina strategies skill not needed (no docking config change)
