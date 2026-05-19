"""build_evomap.py --- Molecular Evolution Map build script.

Reads iteration_log.jsonl + result.log, produces a self-contained
D3.js HTML tree showing Agent hypothesis evolution with RDKit JS
molecule rendering.

Usage:
    python3 tools/build_evomap.py
    python3 tools/build_evomap.py --output output/my_evomap.html
"""

import json
import re
from pathlib import Path
from typing import Any


# ── Data parsing ────────────────────────────────────────────────────────────


def parse_iteration_log(path: Path) -> list[dict[str, Any]]:
    """Parse cumulative iteration_log.jsonl into a list of round entries.

    Each entry: {round, hypothesis_id, success, summary, timestamp, best_be, ...}
    best_be is None at parse time; filled later by merge with result.log data.
    """
    entries: list[dict[str, Any]] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            entries.append({
                "round": obj.get("round"),
                "hypothesis_id": obj.get("hypothesis_id", "UNKNOWN"),
                "success": obj.get("success", False),
                "summary": obj.get("summary", ""),
                "timestamp": obj.get("timestamp", ""),
                "best_be": None,
                "avg_be": None,
                "trivial_count": 0,
                "molecule_count": 0,
                "molecules": [],
            })
    return entries


def parse_result_log(path: Path) -> list[dict[str, Any]]:
    """Parse structured result.log, grouping events by run (delimited by 'start').

    Returns list of runs, each: {round, start_ts, metrics: {...}, molecules: [...]}.
    Non-event entries (stdout, stage) are silently ignored.
    """
    runs: list[dict[str, Any]] = []
    current_run: dict[str, Any] | None = None

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue

            etype = obj.get("type")
            if etype is None:
                continue

            if etype == "start":
                if current_run is not None:
                    runs.append(current_run)
                current_run = {
                    "round": obj.get("round", 0),
                    "start_ts": obj.get("timestamp", ""),
                    "metrics": {},
                    "molecules": [],
                }
            elif etype == "metrics" and current_run is not None:
                trivial = obj.get("trivial_count", 0)
                total = obj.get("molecule_count", 0)
                current_run["metrics"] = {
                    "molecule_count": total,
                    "non_trivial_count": obj.get("non_trivial_count", 0),
                    "trivial_count": trivial,
                    "trivial_ratio": trivial / total if total > 0 else 0,
                    "avg_binding_energy": obj.get("avg_binding_energy"),
                    "min_binding_energy": obj.get("min_binding_energy"),
                }
            elif etype == "molecule" and current_run is not None:
                current_run["molecules"].append({
                    "smiles": obj.get("mol_smiles", ""),
                    "be": obj.get("binding_energy"),
                    "qed": obj.get("qed"),
                    "trivial": obj.get("trivial", False),
                    "syn_steps": obj.get("syn_steps"),
                    "route_quality": obj.get("route_quality"),
                    "composite_score": obj.get("composite_score"),
                })

    if current_run is not None:
        runs.append(current_run)

    return runs


def parse_molecules_jsonl(path: Path) -> list[dict[str, Any]]:
    """Parse output/molecules.jsonl — each line a pipeline run with molecule array.

    Returns list of runs: {timestamp, molecules: [{smiles, be, qed, trivial}, ...]}
    """
    runs: list[dict[str, Any]] = []
    if not path.exists():
        return runs
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                if obj.get("molecules"):
                    runs.append(obj)
            except json.JSONDecodeError:
                continue
    return runs


# ── Tree construction ────────────────────────────────────────────────────────


def _extract_best_be(summary: str) -> float | None:
    """Extract best binding energy from summary text.

    Finds all negative floats in BE-adjacent contexts and returns
    the most negative (best) value.
    """
    # Collect candidates from kcal/mol contexts and all negative floats
    candidates: list[float] = []

    # Floats directly followed by "kcal/mol"
    for v in re.findall(r'([-]?\d+\.\d+)\s*kcal/mol', summary):
        candidates.append(float(v))

    # All negative floats (BE is always negative for Vina)
    if not candidates:
        for v in re.findall(r'-(\d+\.\d+)', summary):
            candidates.append(-float(v))

    if candidates:
        return min(candidates)
    return None


def build_tree(
    iter_entries: list[dict[str, Any]],
    runs: list[dict[str, Any]],
    mol_runs: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Merge iteration log entries with run metrics and molecule data.

    Matching: runs aligned to iter entries by round number + chronological order.
    Molecule data from molecules.jsonl matched by timestamp order to entries.
    Parentage: each entry's parent is the most recent ACCEPTED hypothesis before it.
    """
    run_by_round: dict[int, list[dict[str, Any]]] = {}
    for run in runs:
        rn = run["round"]
        run_by_round.setdefault(rn, []).append(run)

    run_consumed: dict[int, int] = {}

    for entry in iter_entries:
        rn = entry["round"]
        candidates = run_by_round.get(rn, [])
        idx = run_consumed.get(rn, 0)
        if idx < len(candidates):
            run = candidates[idx]
            run_consumed[rn] = idx + 1
            m = run.get("metrics", {})
            entry["best_be"] = m.get("min_binding_energy") or _extract_best_be(entry["summary"])
            entry["avg_be"] = m.get("avg_binding_energy")
            entry["trivial_count"] = m.get("trivial_count", 0)
            entry["molecule_count"] = m.get("molecule_count", 0)
            mols = sorted(run.get("molecules", []),
                         key=lambda x: x.get("be") or 999)
            entry["molecules"] = mols[:5]
        else:
            entry["best_be"] = _extract_best_be(entry["summary"])

    # Fallback: assign unmatched runs to entries lacking molecule data
    all_consumed = sum(run_consumed.values())
    all_available = sum(len(v) for v in run_by_round.values())
    if all_consumed < all_available:
        unmatched = []
        for rn, lst in run_by_round.items():
            used = run_consumed.get(rn, 0)
            unmatched.extend(lst[used:])
        unmatched.sort(key=lambda r: r.get("start_ts", ""))
        entries_need = [e for e in iter_entries if not e["molecules"]]
        entries_need.reverse()
        for i, run in enumerate(unmatched):
            if i >= len(entries_need):
                break
            entry = entries_need[i]
            m = run.get("metrics", {})
            entry["best_be"] = m.get("min_binding_energy") or entry.get("best_be") or _extract_best_be(entry["summary"])
            entry["avg_be"] = m.get("avg_binding_energy") or entry.get("avg_be")
            entry["trivial_count"] = m.get("trivial_count", 0)
            entry["molecule_count"] = m.get("molecule_count", 0)
            mols = sorted(run.get("molecules", []),
                         key=lambda x: x.get("be") or 999)
            entry["molecules"] = mols[:5]

    # Merge molecule data from molecules.jsonl (matched by timestamp order)
    if mol_runs:
        mol_runs_sorted = sorted(mol_runs, key=lambda r: r.get("timestamp", ""))
        entries_need_mols = [e for e in iter_entries if not e["molecules"]]
        entries_need_mols.reverse()
        for i, mr in enumerate(mol_runs_sorted):
            if i >= len(entries_need_mols):
                break
            entry = entries_need_mols[i]
            mols = sorted(mr.get("molecules", []),
                         key=lambda x: x.get("be") or 999)
            entry["molecules"] = mols[:5]
            entry["molecule_count"] = len(mr["molecules"])

    last_accepted: dict[int, str] = {}
    for entry in iter_entries:
        rn = entry["round"]
        parent = None
        for pr in sorted(last_accepted.keys(), reverse=True):
            if pr < rn or (pr == rn and last_accepted[pr] != entry["hypothesis_id"]):
                parent = last_accepted[pr]
                break
        entry["parent"] = parent
        if entry["success"]:
            last_accepted[rn] = entry["hypothesis_id"]

    return iter_entries


# ── HTML generation ──────────────────────────────────────────────────────────


HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Molecular Evolution Map --- MolCraft Agent</title>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#1a1a2e;color:#eee;overflow:hidden;height:100vh}
#tree-container{width:100%;height:100%;position:relative}

.node circle{stroke-width:2px}
.node .success circle{fill:#2ecc71;stroke:#27ae60}
.node .failure circle{fill:#95a5a6;stroke:#7f8c8d}
.node text{font-size:12px;fill:#ecf0f1;font-family:monospace}
.node .be-label{font-size:10px;fill:#bdc3c7}

.link{fill:none;stroke-width:2px}
.link.success{stroke:#27ae60}
.link.failure{stroke:#7f8c8d;stroke-dasharray:6,4}

.tooltip{
  position:absolute;padding:10px 14px;background:rgba(44,62,80,0.95);
  border:1px solid #34495e;border-radius:6px;font-size:13px;
  pointer-events:none;max-width:360px;line-height:1.5;
  box-shadow:0 4px 12px rgba(0,0,0,0.4)
}
.tooltip .hid{font-weight:bold;color:#3498db;font-size:15px}
.tooltip .bebig{color:#e74c3c;font-size:16px;font-weight:bold}

#detail-panel{
  position:absolute;bottom:0;left:0;right:0;height:200px;
  background:rgba(30,39,56,0.97);border-top:1px solid #34495e;
  overflow-x:auto;overflow-y:hidden;white-space:nowrap;
  padding:12px 16px;display:none
}
#detail-panel .mol-card{
  display:inline-block;width:180px;margin-right:12px;
  background:#243447;border-radius:6px;padding:8px;
  text-align:center;vertical-align:top
}
.mol-card .be{color:#e74c3c;font-weight:bold}
.mol-card .qed{color:#bdc3c7;font-size:11px}
.mol-card .route-badge{
  display:inline-block;padding:2px 6px;border-radius:3px;font-size:10px;margin-top:4px
}
.mol-card .route-badge.trivial{background:#7f8c8d;color:#fff}
.mol-card .route-badge.non-trivial{background:#27ae60;color:#fff}
.empty-state{text-align:center;padding:40px;color:#7f8c8d;font-size:16px}
</style>
</head>
<body>
<div id="tree-container">
  <svg width="100%" height="100%"></svg>
  <div class="tooltip" style="display:none"></div>
</div>
<div id="detail-panel"><div style="padding:8px;color:#bdc3c7;font-size:12px">Click a node to see molecule details</div></div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
document.addEventListener("DOMContentLoaded",() => {
  const DATA = __DATA_PLACEHOLDER__;

  if (!DATA || DATA.length === 0) {
    document.getElementById("tree-container").innerHTML =
      '<div class="empty-state">No evolution data yet<br><small>Run main.py first, then rebuild with build_evomap.py</small></div>';
    return;
  }

  const svg = d3.select("#tree-container svg");
  const tooltip = d3.select(".tooltip");
  const detailPanel = d3.select("#detail-panel");
  const margin = {top:40, right:120, bottom:40, left:160};

  const rootData = DATA.find(d => d.parent === null) || DATA[0];
  const stratify = d3.stratify()
    .id(d => d.hypothesis_id)
    .parentId(d => d.parent);
  const root = stratify(DATA);
  root.each(d => { d._children = d.children; });

  const width = window.innerWidth;
  const height = window.innerHeight - 200;
  const treeLayout = d3.tree()
    .size([height - margin.top - margin.bottom, width - margin.left - margin.right]);
  treeLayout(root);

  const g = svg.append("g").attr("transform",`translate(${margin.left},${margin.top})`);

  // Edges
  g.selectAll(".link")
    .data(root.links())
    .join("path")
    .attr("class", d => `link ${d.target.data.success ? "success" : "failure"}`)
    .attr("d", d3.linkHorizontal().x(d => d.y).y(d => d.x))
    .append("title")
    .text(d => `best BE: ${d.target.data.best_be != null ? d.target.data.best_be.toFixed(2) : "N/A"}`);

  // Nodes
  const node = g.selectAll(".node")
    .data(root.descendants())
    .join("g")
    .attr("class", d => `node ${d.data.success ? "success" : "failure"}`)
    .attr("transform", d => `translate(${d.y},${d.x})`)
    .on("mouseover", (event, d) => {
      const dd = d.data;
      tooltip.style("display","block")
        .html(`<div class="hid">${dd.hypothesis_id}</div>
               <div>Best BE: <span class="bebig">${dd.best_be != null ? dd.best_be.toFixed(3) + " kcal/mol" : "N/A"}</span></div>
               <div>Avg BE: ${dd.avg_be != null ? dd.avg_be.toFixed(3) + " kcal/mol" : "N/A"}</div>
               <div>Mols: ${dd.molecule_count || "N/A"} | Trivial: ${dd.trivial_count}</div>
               <div style="margin-top:6px;color:#bdc3c7;font-size:11px;max-height:60px;overflow:hidden">${(dd.summary||"").substring(0,200)}</div>`);
    })
    .on("mousemove", (event) => {
      tooltip.style("left",(event.pageX+12)+"px").style("top",(event.pageY-28)+"px");
    })
    .on("mouseout", () => tooltip.style("display","none"))
    .on("click", (event, d) => showMolecules(d.data));

  node.append("circle")
    .attr("r", d => d.data.hypothesis_id === "BASELINE" ? 8 : 6);

  node.append("text")
    .attr("dy", -12)
    .attr("text-anchor","middle")
    .text(d => d.data.hypothesis_id);

  node.filter(d => d.data.best_be != null).append("text")
    .attr("dy", 20)
    .attr("text-anchor","middle")
    .attr("class","be-label")
    .text(d => d.data.best_be.toFixed(1));

  function showMolecules(data) {
    const mols = data.molecules || [];
    if (mols.length === 0) {
      detailPanel.style("display","block")
        .html(`<div style="padding:12px;color:#7f8c8d">${data.hypothesis_id}: no molecule data saved</div>`);
      return;
    }
    let cards = `<div style="padding:4px 0 4px 12px;color:#bdc3c7;font-size:12px">${data.hypothesis_id} --- ${(data.summary||"").substring(0,120)}</div>`;
    mols.forEach(m => {
      const isTriv = m.trivial;
      cards += `<div class="mol-card">
        <div id="mol-${m.smiles.replace(/[^a-zA-Z0-9]/g,'')}" style="width:160px;height:100px;margin:0 auto"></div>
        <div class="be">BE: ${m.be != null ? m.be.toFixed(2) : "N/A"}</div>
        <div class="qed">QED: ${m.qed != null ? m.qed.toFixed(2) : "N/A"} | Steps: ${m.syn_steps || "?"}</div>
        <span class="route-badge ${isTriv?'trivial':'non-trivial'}">${isTriv?'trivial':'valid route'}</span>
      </div>`;
    });
    detailPanel.style("display","block").html(cards);

    if (typeof initRDKit === "undefined") {
      const s = document.createElement("script");
      s.src = "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js";
      s.onload = () => initRDKit().then(() => renderMols(mols));
      document.head.appendChild(s);
    } else {
      initRDKit().then(() => renderMols(mols));
    }
  }

  function renderMols(mols) {
    mols.forEach(m => {
      const divId = "mol-" + m.smiles.replace(/[^a-zA-Z0-9]/g,'');
      const el = document.getElementById(divId);
      if (!el) return;
      try {
        const mol = RDKitModule.get_mol(m.smiles);
        if (!mol) return;
        mol.draw_to_canvas(el, 160, 100);
        mol.delete();
      } catch(e) {
        el.innerHTML = '<span style="color:#e74c3c;font-size:10px">render error</span>';
      }
    });
  }

  window.addEventListener("resize", () => {
    svg.attr("width", window.innerWidth).attr("height", window.innerHeight);
  });
});
</script>
</body>
</html>
'''


def generate_html(tree: list[dict[str, Any]], output_path: Path) -> None:
    """Embed tree data as JSON into the HTML template and write to output_path."""
    data_json = json.dumps(
        [{
            "hypothesis_id": node["hypothesis_id"],
            "round": node["round"],
            "success": node["success"],
            "parent": node["parent"],
            "best_be": node.get("best_be"),
            "avg_be": node.get("avg_be"),
            "trivial_count": node.get("trivial_count", 0),
            "molecule_count": node.get("molecule_count", 0),
            "summary": node.get("summary", ""),
            "molecules": node.get("molecules", []),
        } for node in tree],
        ensure_ascii=False,
        allow_nan=False,
    )
    html = HTML_TEMPLATE.replace("__DATA_PLACEHOLDER__", data_json)
    output_path.write_text(html, encoding="utf-8")


# ── CLI ──────────────────────────────────────────────────────────────────────


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Build Molecular Evolution Map HTML"
    )
    parser.add_argument(
        "--iteration-log", type=Path, default=Path("docs/iteration_log.jsonl"),
        help="Path to iteration_log.jsonl",
    )
    parser.add_argument(
        "--result-log", type=Path, default=Path("output/result.log"),
        help="Path to result.log",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("output/evomap.html"),
        help="Output HTML path",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    iter_log = args.iteration_log if args.iteration_log.is_absolute() else project_root / args.iteration_log
    result_log = args.result_log if args.result_log.is_absolute() else project_root / args.result_log
    output_path = args.output if args.output.is_absolute() else project_root / args.output

    if not iter_log.exists():
        print(f"iteration_log not found: {iter_log}")
        raise SystemExit(1)

    print(f"Reading iteration log: {iter_log}")
    iter_entries = parse_iteration_log(iter_log)
    print(f"  {len(iter_entries)} entries")

    runs: list[dict[str, Any]] = []
    if result_log.exists():
        print(f"Reading result log: {result_log}")
        runs = parse_result_log(result_log)
        print(f"  {len(runs)} pipeline runs")
    else:
        print("result.log not found, tree nodes will lack molecule details")

    print("Building tree...")
    mol_runs = parse_molecules_jsonl(project_root / "output" / "molecules.jsonl")
    if mol_runs:
        print(f"  {len(mol_runs)} molecule records loaded")
    tree = build_tree(iter_entries, runs, mol_runs)

    print(f"Generating HTML: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    generate_html(tree, output_path)

    n_success = sum(1 for n in tree if n["success"])
    n_fail = sum(1 for n in tree if not n["success"])
    print(f"\nDone! Open in browser: {output_path}")
    print(f"  Nodes: {len(tree)}")
    print(f"  Accepted: {n_success} | Rejected: {n_fail}")


if __name__ == "__main__":
    main()
