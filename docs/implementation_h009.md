# H009 — Consensus Docking Implementation

## Changes Made

### 1. `src/docking.py`

#### `dock_molecule()` — added `seed` parameter
- Now accepts optional `seed` parameter (int, default 0 → uses 42).
- Backward compatible: all existing callers (`batch_dock`, `pipeline.py`) continue to work unchanged.
- Internal: `vina_seed = seed if seed != 0 else 42`

#### `dock_molecule_consensus()` — new function
- Runs N independent dockings (default 3) with different seeds [42, 123, 456].
- Returns median `binding_energy`, `std_energy` (stdev), `all_energies` list.
- Graceful degradation: if only 1/3 succeeds, returns that single result.
- Requires ≥2 successful dockings for a consensus result.

### 2. `tools/pipeline.py`

#### Added import
- `from docking import dock_molecule_consensus`

#### Consensus re-ranking block (after dedup, before final selection)
- Takes top `n_top * 2` candidates (max 20) after evolutionary generations.
- Re-scores each with `dock_molecule_consensus()`.
- Adds `consensus_std` and `consensus_n` fields to result dicts.
- Replaces single-run `binding_energy` with median consensus energy.
- Re-sorts by consensus score, selects final top-N.
- Logs per-molecule consensus results with std for transparency.

## Design Decisions

| Decision | Rationale |
|----------|-----------|
| Consensus only on final top-20, not during evolution | Cost: 20×2 extra dockings (~7 min) vs 50×3×2 (~50 min) for full consensus |
| Median (not mean) | Robust to outliers from failed docking poses |
| 3 runs with fixed seeds [42,123,456] | Reproducible; 3 is the minimum for meaningful median |
| Fallback to single-run on partial failure | Prevents losing good candidates if 1/3 fails |
