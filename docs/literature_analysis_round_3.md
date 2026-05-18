# Literature Analysis — Round 3

## Re-read: Coscientist (Boiko et al., 2023)

### Key Patterns Extracted

| # | Pattern | Coscientist Implementation | Mappable to MolCraft |
|---|---------|---------------------------|-----------------------|
| 1 | **Self-error-correction** | Agent detects failed experiments and re-runs with adjusted parameters | Already partially done (SA filter, retrosynthesis retry) |
| 2 | **Consensus scoring** | Performs experiments multiple times, uses average/consensus result | **H009: Consensus docking** ← selected |
| 3 | **Web-search iteration** | Searches literature iteratively to refine protocols | Not directly applicable (no web interface) |
| 4 | **Modular code generation** | Generates and executes Python code blocks | Partially done (code is pre-built) |

### Why Consensus Docking?

- **Directly mappable**: Vina docking is stochastic — different random seeds produce different energy scores (±0.1–0.5 kcal/mol for the same molecule).
- **Low risk**: The consensus function wraps `dock_molecule()`, no changes to the core algorithm.
- **Addresses the scoring reliability problem**: In Round 2, top-10 rankings were dominated by sulfonamides partly because single-run docking noise favored certain scaffolds.
- **Cost**: 3 runs × top-20 candidates = ~60 extra dockings ≈ 10–15 minutes extra runtime. Acceptable in the final (3rd) round.

### Other Papers Considered

- **MolLEO (Wang et al., 2024b)**: Multi-objective optimization with LLM operators. Already implemented (evolutionary pipeline with mutation).
- **MOOSE-Chem (Yang et al., 2025)**: Evolutionary algorithms for chemical space. Already the backbone of our pipeline.
- **Deep Docking (Gentile et al., 2020)**: Virtual screening acceleration. Not applicable (we generate de novo, not screen libraries).

## Selected Hypothesis: H009 — Consensus Docking
