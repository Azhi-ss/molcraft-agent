# Diagnosis — Round 3

## Problem Statement

After Rounds 1–2, the pipeline produces chemically diverse molecules but the top-10 ranking is dominated by sulfonamide scaffolds. The hypothesis is that *single-run Vina scores contain stochastic noise* that biases ranking toward certain scaffolds.

## Root Cause: Stochastic Docking Noise

- Vina uses a stochastic Monte Carlo search algorithm (seed=42 hardcoded).
- Different seeds produce slightly different binding energies for identical molecules (±0.1–0.5 kcal/mol).
- When 50+ molecules are compared on single-run scores, the 0.1–0.5 kcal/mol noise can shuffle the ranking of truly similar binders.
- This biases selection toward scaffolds that happen to score well in a single run, reducing chemical diversity in the final set.

## Solution: H009 — Consensus Docking

### Implementation
- `dock_molecule_consensus()` runs 3 independent dockings with seeds (42, 123, 456).
- Takes the median energy as consensus score.
- Applied to the top-20 candidates after evolutionary generations (not during them — too expensive).
- Evolutions still use fast single-run docking; only final selection uses consensus.

### Expected Impact
- More reliable final ranking → genuinely better binders selected.
- Reduced scaffold bias → potentially more diverse top-10.
- Minimal code risk: pure additive change.

## Metrics to Watch
- **Best binding energy**: should improve or stay same
- **Average binding energy**: should improve or stay same
- **Chemical diversity**: broader scaffold representation in top-10
- **Consensus std**: indicates scoring reliability
