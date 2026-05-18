# Molecular Evolution Map — Design Spec

> **Goal:** Visualize Agent's research iteration path across rounds, showing hypothesis lineage and molecule outputs in a single interactive HTML.

**Architecture:** Offline build script reads existing JSONL logs, produces a self-contained D3.js HTML tree. Zero pipeline changes, zero runtime dependencies beyond a browser.

**Tech Stack:** Python 3 (build), D3.js v7 CDN (tree), RDKit JS CDN (molecule 2D rendering)

---

## Architecture

```
docs/iteration_log.jsonl --+
                            +--> tools/build_evomap.py --> output/evomap.html
output/result.log ---------+
```

Three units, each single-file:

| Unit | Path | Role |
|------|------|------|
| Build script | `tools/build_evomap.py` | Read both JSONL sources, merge into one JSON tree, embed into HTML template |
| HTML template | inline in build script | Static HTML shell with `<script>` that loads data and renders tree |
| Output | `output/evomap.html` | Self-contained deliverable, open in browser |

## Data Model

### Input: iteration_log.jsonl

Each line is a manually-reported iteration record:
```json
{"round": 7, "hypothesis_id": "H015", "success": true, "summary": "..."}
```

### Input: result.log

Structured JSONL with typed events. Relevant events:
- `type: "start"` — marks a new pipeline run, contains `round`
- `type: "molecule"` — individual molecule: `mol_smiles`, `binding_energy`, `composite_score`, `syn_steps`, `qed`, `trivial`, `route_quality`
- `type: "metrics"` — aggregate: `molecule_count`, `non_trivial_count`, `trivial_count`, `avg_binding_energy`, `min_binding_energy`
- `type: "stage"` — marks pipeline phases

### Output tree node (merged)

```json
{
  "id": "BASELINE",
  "round": 1,
  "success": true,
  "parent": null,
  "best_be": -8.56,
  "avg_be": -8.058,
  "trivial_ratio": 0.0,
  "summary": "...",
  "molecules": [
    {"smiles": "Cc1ccc...", "be": -8.56, "qed": 0.72, "trivial": false}
  ]
}
```

### Parentage inference

Nodes link by round order: round N's parent is round N-1. Hypothesis lineage is reconstructed from the chronological sequence in iteration_log.jsonl. REJECTED hypotheses (H004, H012) branch off and dead-end — marked with dashed edges and gray nodes.

## Tree Visualization

### Layout

- D3.js tree layout, horizontal (root at left, branches expand right)
- Root node: "BASELINE" with the initial run's best BE
- Each child: one hypothesis (H009, H010, H011, ...)
- Leaf nodes: most recent accepted hypotheses

### Visual encoding

| Element | Success | Failure (REJECTED) |
|---------|---------|---------------------|
| Node fill | Green (#2ecc71) | Gray (#bdc3c7) |
| Edge style | Solid, dark | Dashed, light gray |
| Edge label | best BE (e.g. "-9.94") | best BE |
| Node label | Hypothesis ID | Hypothesis ID |

### Interaction

- **Hover node**: tooltip with hypothesis ID, best BE, avg BE, trivial ratio, round summary
- **Click node**: expands a detail panel below the tree showing:
  - Round summary text
  - Top molecules rendered as 2D structures (RDKit JS `MolDraw2D`)
  - Each molecule card: structure image, BE, QED, synthesis steps

## Build Script Logic

```
1. Parse iteration_log.jsonl → list of {round, hypothesis_id, success, summary}
2. Parse result.log → group MoleculeEvents by run (using "start" events as delimiters)
3. Match runs to hypothesis rounds by chronological alignment
4. Build tree: root = earliest baseline run, children = subsequent experiments
5. Generate HTML with JSON data embedded as `const DATA = {...};`
6. Write output/evomap.html
```

## Scope boundaries

**Included:**
- Tree of hypothesis rounds with BE and trivial metrics
- Molecule 2D structure rendering on click
- Success/failure visual distinction
- Self-contained HTML output

**Excluded:**
- SmartSMILES — using standard SMILES + RDKit JS rendering
- Pipeline modifications — purely offline processing
- Multi-run lineage tracking within a single generation (generation-level parentage)
- Server-side rendering or WebSocket updates
- Mobile responsive layout (desktop-first, tree needs width)
