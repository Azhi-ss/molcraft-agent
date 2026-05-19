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

    # Unique IDs: d3.stratify requires unique ids.
    id_counts: dict[str, int] = {}
    for e in iter_entries:
        hid = e["hypothesis_id"]
        id_counts[hid] = id_counts.get(hid, 0) + 1
    id_counter: dict[str, int] = {}
    for e in iter_entries:
        hid = e["hypothesis_id"]
        if id_counts[hid] > 1:
            id_counter[hid] = id_counter.get(hid, 0) + 1
            new_id = f"{hid}.{id_counter[hid]}"
            # Update parent refs that pointed to this (original) id
            for other in iter_entries:
                if other.get("parent") == hid:
                    other["parent"] = new_id
            e["hypothesis_id"] = new_id

    # Ensure single root for D3 tree: nodes with parent=null → virtual root
    root_count = sum(1 for e in iter_entries if e["parent"] is None)
    if root_count > 1:
        virtual_root = {
            "round": 0,
            "hypothesis_id": "MOLCRAFT",
            "success": True,
            "summary": "MolCraft Agent",
            "timestamp": "",
            "best_be": None,
            "avg_be": None,
            "trivial_count": 0,
            "molecule_count": 0,
            "parent": None,
            "molecules": [],
        }
        for e in iter_entries:
            if e["parent"] is None:
                e["parent"] = "MOLCRAFT"
        iter_entries.insert(0, virtual_root)

    return iter_entries


# -- HTML generation --

# Cytoscape.js + comparison table template.
# Two placeholders: __GRAPH_PLACEHOLDER__ (Cytoscape elements JSON)
# and __SESSIONS_PLACEHOLDER__ (session comparison table data)

HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Molecular Evolution Map -- MolCraft Agent</title>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#1a1a2e;color:#eee;height:100vh;display:flex;flex-direction:column}
#cy{flex:0 0 60vh;background:#1a1a2e;border-bottom:2px solid #34495e}
#lower{flex:1;display:flex;flex-direction:column;overflow:hidden}
#tabs{display:flex;background:#16213e;border-bottom:1px solid #34495e}
.tab{padding:8px 18px;cursor:pointer;font-size:13px;color:#7f8c8d;border-bottom:2px solid transparent}
.tab.active{color:#3498db;border-bottom-color:#3498db}
#table-wrap{flex:1;overflow:auto;padding:8px 12px}
#mol-panel{flex:1;overflow-x:auto;overflow-y:hidden;white-space:nowrap;padding:12px 16px;display:none}
#mol-panel .mol-card{display:inline-block;width:180px;margin-right:12px;background:#243447;border-radius:6px;padding:8px;text-align:center;vertical-align:top}
.mol-card .be{color:#e74c3c;font-weight:bold}
.mol-card .qed{color:#bdc3c7;font-size:11px}
.mol-card .route-badge{display:inline-block;padding:2px 6px;border-radius:3px;font-size:10px;margin-top:4px}
.mol-card .route-badge.trivial{background:#7f8c8d;color:#fff}
.mol-card .route-badge.non-trivial{background:#27ae60;color:#fff}

table{border-collapse:collapse;width:100%;font-size:12px}
th,td{padding:6px 10px;border:1px solid #2c3e50;text-align:center;white-space:nowrap}
th{background:#16213e;position:sticky;top:0;z-index:1}
tr:nth-child(even){background:rgba(44,62,80,0.3)}
.session-col{color:#7f8c8d;font-size:10px}
.be-up{color:#2ecc71}
.be-down{color:#e74c3c}
.empty-state{text-align:center;padding:40px;color:#7f8c8d;font-size:16px}

.tooltip{position:absolute;padding:10px 14px;background:rgba(44,62,80,0.95);border:1px solid #34495e;border-radius:6px;font-size:13px;pointer-events:none;max-width:360px;line-height:1.5;box-shadow:0 4px 12px rgba(0,0,0,0.4);display:none;z-index:10}
.tooltip .hid{font-weight:bold;color:#3498db;font-size:15px}
.tooltip .bebig{color:#e74c3c;font-size:16px;font-weight:bold}

#legend{position:absolute;top:8px;right:12px;background:rgba(22,33,62,0.9);padding:8px 12px;border-radius:6px;font-size:11px;z-index:10}
.legend-dot{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:4px;vertical-align:middle}
</style>
</head>
<body>
<div id="cy"></div>
<div id="lower">
  <div id="tabs">
    <div class="tab active" data-tab="table">Session Compare</div>
    <div class="tab" data-tab="mol">Molecules</div>
  </div>
  <div id="table-wrap"><table id="session-table"></table></div>
  <div id="mol-panel"><div style="padding:8px;color:#bdc3c7;font-size:12px">Click a graph node to see molecule structures</div></div>
</div>
<div id="legend"></div>
<div class="tooltip"></div>

<script src="cytoscape.min.js"></script>
<script>
(function() {
  var GRAPH = __GRAPH_PLACEHOLDER__;
  var SESSIONS = __SESSIONS_PLACEHOLDER__;

  if (!GRAPH || GRAPH.length === 0) {
    document.getElementById("cy").innerHTML =
      '<div class="empty-state">No data yet<br><small>Run main.py, then rebuild with build_evomap.py</small></div>';
    return;
  }

  var sessionColors = ["#3498db","#e67e22","#2ecc71","#9b59b6","#1abc9c","#f39c12"];
  var sessionColorMap = {};
  SESSIONS.forEach(function(s, i) {
    sessionColorMap[s] = sessionColors[i % sessionColors.length];
  });

  var cy = cytoscape({
    container: document.getElementById("cy"),
    elements: GRAPH,
    style: [
      { selector: "node", style: {
        "label": "data(id)",
        "color": "#ecf0f1",
        "font-size": "11px",
        "font-family": "monospace",
        "text-valign": "center",
        "text-halign": "center",
        "background-color": "data(bg)",
        "border-color": "data(border)",
        "border-width": 2.5,
        "width": "data(size)",
        "height": "data(size)"
      }},
      { selector: "edge", style: {
        "width": 2,
        "line-color": "data(ecolor)",
        "target-arrow-color": "data(ecolor)",
        "target-arrow-shape": "triangle",
        "curve-style": "bezier",
        "line-style": "data(estyle)"
      }}
    ],
    layout: { name: "breadthfirst", directed: true, spacingFactor: 1.3, avoidOverlap: true },
    wheelSensitivity: 0.3,
    maxZoom: 3,
    minZoom: 0.3
  });

  // Tooltip
  var tooltip = document.querySelector(".tooltip");
  cy.on("mouseover", "node", function(evt) {
    var d = evt.target.data();
    tooltip.style.display = "block";
    tooltip.innerHTML =
      '<div class="hid">' + d.id + '</div>' +
      '<div>Best BE: <span class="bebig">' + (d.best_be != null ? d.best_be.toFixed(3) + " kcal/mol" : "N/A") + '</span></div>' +
      '<div>Avg BE: ' + (d.avg_be != null ? d.avg_be.toFixed(3) + " kcal/mol" : "N/A") + '</div>' +
      '<div>Session: ' + (d.run_id || "N/A") + '</div>' +
      '<div style="margin-top:6px;color:#bdc3c7;font-size:11px;max-height:60px;overflow:hidden">' + (d.summary||"").substring(0,200) + '</div>';
    tooltip.style.left = (evt.originalEvent.pageX + 12) + "px";
    tooltip.style.top = (evt.originalEvent.pageY - 28) + "px";
  });
  cy.on("mouseout", "node", function() { tooltip.style.display = "none"; });

  // Click node -> show molecules
  cy.on("click", "node", function(evt) {
    var d = evt.target.data();
    showMolecules(d);
    document.querySelector(".tab[data-tab='mol']").click();
  });

  function showMolecules(data) {
    var mols = data.molecules || [];
    var panel = document.getElementById("mol-panel");
    if (mols.length === 0) {
      panel.innerHTML = '<div style="padding:12px;color:#7f8c8d">' + data.id + ': no molecule data</div>';
      return;
    }
    var cards = '<div style="padding:4px 0;color:#bdc3c7;font-size:12px">' + data.id + '</div>';
    mols.forEach(function(m) {
      var isTriv = m.trivial;
      cards += '<div class="mol-card">' +
        '<div id="mol-' + m.smiles.replace(/[^a-zA-Z0-9]/g,'') + '" style="width:160px;height:100px;margin:0 auto"></div>' +
        '<div class="be">BE: ' + (m.be != null ? m.be.toFixed(2) : "N/A") + '</div>' +
        '<div class="qed">QED: ' + (m.qed != null ? m.qed.toFixed(2) : "N/A") + ' | Steps: ' + (m.syn_steps || "?") + '</div>' +
        '<span class="route-badge ' + (isTriv?"trivial":"non-trivial") + '">' + (isTriv?"trivial":"valid route") + '</span>' +
        '</div>';
    });
    panel.innerHTML = cards;

    if (typeof initRDKit === "undefined") {
      var s = document.createElement("script");
      s.src = "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js";
      s.onload = function() { initRDKit().then(function() { renderMols(mols); }); };
      document.head.appendChild(s);
    } else {
      initRDKit().then(function() { renderMols(mols); });
    }
  }
  function renderMols(mols) {
    mols.forEach(function(m) {
      var divId = "mol-" + m.smiles.replace(/[^a-zA-Z0-9]/g,"");
      var el = document.getElementById(divId);
      if (!el) return;
      try {
        var mol = RDKitModule.get_mol(m.smiles);
        if (!mol) return;
        mol.draw_to_canvas(el, 160, 100);
        mol.delete();
      } catch(e) {
        el.innerHTML = '<span style="color:#e74c3c;font-size:10px">render error</span>';
      }
    });
  }

  // Tab switching
  document.querySelectorAll(".tab").forEach(function(t) {
    t.addEventListener("click", function() {
      document.querySelectorAll(".tab").forEach(function(x) { x.classList.remove("active"); });
      t.classList.add("active");
      var tab = t.dataset.tab;
      document.getElementById("table-wrap").style.display = tab === "table" ? "block" : "none";
      document.getElementById("mol-panel").style.display = tab === "mol" ? "block" : "none";
    });
  });

  // Session comparison table
  buildTable();

  function buildTable() {
    var hypoMap = {};
    GRAPH.forEach(function(el) {
      if (el.group !== "nodes") return;
      var id = el.data.id;
      if (id === "MOLCRAFT") return;
      if (!hypoMap[id]) hypoMap[id] = {};
      hypoMap[id][el.data.run_id || "?"] = el.data;
    });

    var hypoIds = Object.keys(hypoMap).sort();
    if (hypoIds.length === 0) return;

    var thead = "<tr><th>Hypothesis</th>";
    SESSIONS.forEach(function(s) {
      thead += '<th class="session-col">' + s.substring(5, 17) + '</th>';
    });
    thead += "</tr>";

    var tbody = "";
    hypoIds.forEach(function(hid) {
      tbody += "<tr><td style='text-align:left;font-family:monospace'>" + hid + "</td>";
      var prevBe = null;
      SESSIONS.forEach(function(s) {
        var d = hypoMap[hid][s];
        if (d && d.best_be != null) {
          var be = d.best_be;
          var cls = "";
          if (prevBe != null) {
            cls = be < prevBe ? "be-up" : (be > prevBe ? "be-down" : "");
          }
          var icon = d.success ? "&#x2714;" : "&#x2718;";
          tbody += "<td class='" + cls + "'>" + icon + " " + be.toFixed(2) + "</td>";
          prevBe = be;
        } else {
          tbody += "<td style='color:#555'>-</td>";
          prevBe = null;
        }
      });
      tbody += "</tr>";
    });

    document.getElementById("session-table").innerHTML = thead + tbody;
  }

  // Legend
  var legendHtml = "";
  SESSIONS.forEach(function(s, i) {
    legendHtml += '<div><span class="legend-dot" style="background:' + sessionColors[i % sessionColors.length] + '"></span> ' + s.substring(5,17) + '</div>';
  });
  legendHtml += '<div style="margin-top:4px"><span class="legend-dot" style="background:#2ecc71"></span> Accepted</div>';
  legendHtml += '<div><span class="legend-dot" style="background:#95a5a6"></span> Rejected</div>';
  document.getElementById("legend").innerHTML = legendHtml;

})();
</script>
</body>
</html>
'''


def generate_html(tree: list[dict[str, Any]], output_path: Path) -> None:
    """Convert tree data to Cytoscape graph elements + session table, write HTML."""
    if not tree:
        html = HTML_TEMPLATE.replace("__GRAPH_PLACEHOLDER__", "[]").replace("__SESSIONS_PLACEHOLDER__", "[]")
        output_path.write_text(html, encoding="utf-8")
        _copy_cytoscape(output_path)
        return

    elements: list[dict] = []
    sessions_set: set[str] = set()
    for n in tree:
        rid = n.get("run_id", "")
        if rid:
            sessions_set.add(rid)

    sessions = sorted(sessions_set, reverse=True)
    session_colors = ["#3498db", "#e67e22", "#2ecc71", "#9b59b6", "#1abc9c", "#f39c12"]
    session_color_map = {s: session_colors[i % len(session_colors)] for i, s in enumerate(sessions)}

    for n in tree:
        hid = n["hypothesis_id"]
        success = n.get("success", False)
        run_id = n.get("run_id", "unknown")
        border = session_color_map.get(run_id, "#34495e")
        bg = "#2ecc71" if success else "#95a5a6"

        node = {
            "data": {
                "id": hid,
                "bg": bg,
                "border": border,
                "size": 14 if hid == "MOLCRAFT" else 12,
                "best_be": n.get("best_be"),
                "avg_be": n.get("avg_be"),
                "summary": n.get("summary", ""),
                "success": success,
                "run_id": run_id,
                "molecules": n.get("molecules", []),
            },
        }
        elements.append(node)

        parent = n.get("parent")
        if parent:
            edge = {
                "data": {
                    "id": parent + "_to_" + hid,
                    "source": parent,
                    "target": hid,
                    "ecolor": "#27ae60" if success else "#7f8c8d",
                    "estyle": "solid" if success else "dashed",
                }
            }
            elements.append(edge)

    graph_json = json.dumps(elements, ensure_ascii=False, allow_nan=False)
    sessions_json = json.dumps(list(sessions), ensure_ascii=False) if sessions else "[]"
    html = HTML_TEMPLATE.replace("__GRAPH_PLACEHOLDER__", graph_json).replace("__SESSIONS_PLACEHOLDER__", sessions_json)
    output_path.write_text(html, encoding="utf-8")
    _copy_cytoscape(output_path)


def _copy_cytoscape(output_path: Path) -> None:
    """Copy Cytoscape.js to output dir for local loading."""
    src = Path(__file__).resolve().parent.parent / "output" / "cytoscape.min.js"
    dst = output_path.parent / "cytoscape.min.js"
    if src.exists() and src != dst:
        import shutil
        shutil.copy2(src, dst)




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
