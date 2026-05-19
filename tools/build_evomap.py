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

# Modern Dashboard Template with Fira fonts and Sky/Slate color scheme.
HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>EvoMap Pro | Molecular Evolution Explorer</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fira+Code:wght@400;600&family=Fira+Sans:wght@300;400;500;600;700&display=swap">
<style>
:root {
  --bg: #020617;
  --card-bg: #0f172a;
  --panel-bg: #1e293b;
  --border: #334155;
  --text: #f8fafc;
  --text-muted: #94a3b8;
  --primary: #38bdf8;
  --success: #22c55e;
  --danger: #ef4444;
  --accent: #818cf8;
  --header-h: 56px;
}

* { margin: 0; padding: 0; box-sizing: border-box; }
body {
  font-family: 'Fira Sans', sans-serif;
  background: var(--bg);
  color: var(--text);
  height: 100vh;
  display: flex;
  flex-direction: column;
  overflow: hidden;
}

/* Header */
header {
  height: var(--header-h);
  background: var(--card-bg);
  border-bottom: 1px solid var(--border);
  display: flex;
  align-items: center;
  padding: 0 20px;
  justify-content: space-between;
  z-index: 100;
  box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1);
}
.logo { display: flex; align-items: center; gap: 10px; font-weight: 700; font-size: 18px; color: var(--primary); }
.stats { display: flex; gap: 24px; font-size: 13px; color: var(--text-muted); }
.stat-item b { color: var(--text); margin-right: 4px; }

/* Main Layout */
main {
  flex: 1;
  display: flex;
  flex-direction: column;
  overflow: hidden;
}

#cy-wrapper {
  flex: 0 0 55vh;
  position: relative;
  background: radial-gradient(circle at center, #0f172a 0%, #020617 100%);
  border-bottom: 1px solid var(--border);
}
#cy { width: 100%; height: 100%; }

#bottom-panel {
  flex: 1;
  display: flex;
  flex-direction: column;
  background: var(--bg);
  overflow: hidden;
}

/* Tabs */
.tabs-header {
  display: flex;
  background: var(--card-bg);
  border-bottom: 1px solid var(--border);
  padding: 0 10px;
}
.tab-btn {
  padding: 14px 24px;
  cursor: pointer;
  font-size: 14px;
  font-weight: 500;
  color: var(--text-muted);
  border-bottom: 2px solid transparent;
  transition: all 0.2s ease;
}
.tab-btn:hover { color: var(--text); background: rgba(255,255,255,0.03); }
.tab-btn.active { color: var(--primary); border-bottom-color: var(--primary); }

.tab-content { flex: 1; overflow: auto; display: none; padding: 20px; }
.tab-content.active { display: block; }

/* Table Styling */
table { width: 100%; border-collapse: separate; border-spacing: 0; font-size: 13px; }
th {
  background: var(--panel-bg);
  padding: 12px 16px;
  text-align: left;
  font-weight: 600;
  color: var(--text-muted);
  position: sticky;
  top: 0;
  border-bottom: 1px solid var(--border);
}
td { padding: 10px 16px; border-bottom: 1px solid rgba(51, 65, 85, 0.5); }
tr:hover { background: rgba(56, 189, 248, 0.05); }
.be-val { font-family: 'Fira Code', monospace; }
.diff-up { color: var(--success); }
.diff-down { color: var(--danger); }

/* Molecule Cards */
.mol-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(220px, 1fr));
  gap: 16px;
}
.mol-card {
  background: var(--card-bg);
  border: 1px solid var(--border);
  border-radius: 12px;
  padding: 16px;
  transition: transform 0.2s ease, border-color 0.2s ease;
  position: relative;
  overflow: hidden;
}
.mol-card:hover { transform: translateY(-4px); border-color: var(--primary); }
.mol-canvas { width: 100%; height: 140px; margin-bottom: 12px; background: white; border-radius: 8px; }
.mol-info { display: flex; flex-direction: column; gap: 4px; }
.mol-be { font-size: 18px; font-weight: 700; color: var(--primary); font-family: 'Fira Code', monospace; }
.mol-meta { font-size: 11px; color: var(--text-muted); display: flex; justify-content: space-between; }
.badge {
  padding: 2px 8px;
  border-radius: 99px;
  font-size: 10px;
  font-weight: 600;
  text-transform: uppercase;
}
.badge-success { background: rgba(34, 197, 94, 0.2); color: #4ade80; }
.badge-muted { background: rgba(148, 163, 184, 0.2); color: #cbd5e1; }

/* Tooltip & Legends */
.tooltip {
  position: absolute;
  padding: 12px 16px;
  background: rgba(15, 23, 42, 0.95);
  backdrop-filter: blur(8px);
  border: 1px solid var(--border);
  border-radius: 10px;
  box-shadow: 0 10px 15px -3px rgb(0 0 0 / 0.5);
  pointer-events: none;
  z-index: 1000;
  display: none;
  max-width: 320px;
}
.tooltip h4 { color: var(--primary); margin-bottom: 4px; font-size: 14px; }
.tooltip .be { font-size: 16px; font-weight: 700; color: var(--text); }

#legend {
  position: absolute;
  bottom: 20px;
  left: 20px;
  background: rgba(15, 23, 42, 0.8);
  backdrop-filter: blur(4px);
  padding: 12px;
  border-radius: 8px;
  border: 1px solid var(--border);
  font-size: 11px;
  z-index: 10;
}
.legend-item { display: flex; align-items: center; gap: 8px; margin-bottom: 4px; }
.dot { width: 10px; height: 10px; border-radius: 50%; }

/* Scrollbar */
::-webkit-scrollbar { width: 8px; height: 8px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb { background: var(--border); border-radius: 4px; }
::-webkit-scrollbar-thumb:hover { background: var(--text-muted); }

.empty-state {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  height: 100%;
  color: var(--text-muted);
  text-align: center;
}
</style>
</head>
<body>

<header>
  <div class="logo">
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"></path><polyline points="3.27 6.96 12 12.01 20.73 6.96"></polyline><line x1="12" y1="22.08" x2="12" y2="12"></line></svg>
    EvoMap Pro
  </div>
  <div class="stats" id="header-stats">
    <!-- Filled by JS -->
  </div>
</header>

<main>
  <div id="cy-wrapper">
    <div id="cy"></div>
    <div id="legend"></div>
  </div>

  <div id="bottom-panel">
    <div class="tabs-header">
      <div class="tab-btn active" data-tab="compare">Evolution History</div>
      <div class="tab-btn" data-tab="molecules">Molecule Gallery</div>
    </div>
    
    <div id="compare" class="tab-content active">
      <table id="session-table"></table>
    </div>
    
    <div id="molecules" class="tab-content">
      <div id="mol-detail-header" style="margin-bottom:16px; font-weight:600; color:var(--primary)">
        Select a node to view molecules
      </div>
      <div class="mol-grid" id="mol-grid">
        <!-- Filled by JS -->
      </div>
    </div>
  </div>
</main>

<div class="tooltip" id="main-tooltip"></div>

<script src="cytoscape.min.js"></script>
<script>
(function() {
  const GRAPH = __GRAPH_PLACEHOLDER__;
  const SESSIONS = __SESSIONS_PLACEHOLDER__;

  if (!GRAPH || GRAPH.length === 0) {
    document.querySelector("main").innerHTML = `
      <div class="empty-state">
        <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" style="margin-bottom:16px; opacity:0.5">
          <circle cx="12" cy="12" r="10"></circle><line x1="12" y1="8" x2="12" y2="12"></line><line x1="12" y1="16" x2="12.01" y2="16"></line>
        </svg>
        <h3>No Data Available</h3>
        <p>Run your evolution pipeline first, then rebuild this map.</p>
      </div>`;
    return;
  }

  // Header Stats
  const nodesCount = GRAPH.filter(e => e.group === 'nodes').length;
  const acceptedCount = GRAPH.filter(e => e.group === 'nodes' && e.data.success).length;
  document.getElementById('header-stats').innerHTML = `
    <div class="stat-item"><b>${nodesCount}</b> Nodes</div>
    <div class="stat-item"><b>${acceptedCount}</b> Accepted</div>
    <div class="stat-item"><b>${SESSIONS.length}</b> Iterations</div>
  `;

  const cy = window.cy = cytoscape({
    container: document.getElementById("cy"),
    elements: GRAPH,
    style: [
      { selector: "node", style: {
        "label": "data(id)",
        "color": "#fff",
        "font-size": "10px",
        "font-family": "Fira Code, monospace",
        "text-valign": "center",
        "text-halign": "center",
        "background-color": "data(bg)",
        "border-color": "data(border)",
        "border-width": 2,
        "width": "data(size)",
        "height": "data(size)",
        "transition-property": "background-color, border-color, width, height",
        "transition-duration": "0.2s"
      }},
      { selector: "node:selected", style: {
        "border-color": "#38bdf8",
        "border-width": 4,
        "width": 36,
        "height": 36
      }},
      { selector: "edge", style: {
        "width": 2,
        "line-color": "data(ecolor)",
        "target-arrow-color": "data(ecolor)",
        "target-arrow-shape": "triangle",
        "curve-style": "bezier",
        "line-style": "data(estyle)",
        "opacity": 0.6
      }}
    ],
    layout: { 
      name: "breadthfirst", 
      directed: true, 
      spacingFactor: 1.1, 
      avoidOverlap: true,
      padding: 50
    },
    wheelSensitivity: 0.2,
    maxZoom: 2,
    minZoom: 0.2
  });

  // Tooltip Logic
  const tooltip = document.getElementById("main-tooltip");
  cy.on("mouseover", "node", (evt) => {
    const d = evt.target.data();
    tooltip.style.display = "block";
    tooltip.innerHTML = `
      <h4>${d.id}</h4>
      <div class="be">Best BE: ${d.best_be != null ? d.best_be.toFixed(3) : 'N/A'}</div>
      <div style="font-size:11px; color:var(--text-muted); margin: 4px 0 8px">
        Session: ${d.run_id}<br>
        Status: ${d.success ? 'Accepted' : 'Rejected'}
      </div>
      <div style="font-size:11px; max-height:80px; overflow:hidden; border-top: 1px solid var(--border); padding-top:8px">
        ${d.summary || 'No summary available.'}
      </div>
    `;
  });
  
  cy.on("mousemove", (evt) => {
    if (tooltip.style.display === "block") {
      tooltip.style.left = (evt.renderedPosition.x + 20) + "px";
      tooltip.style.top = (evt.renderedPosition.y - 20) + "px";
    }
  });

  cy.on("mouseout", "node", () => { tooltip.style.display = "none"; });

  // Click Interaction
  cy.on("click", "node", (evt) => {
    const d = evt.target.data();
    updateMoleculeGallery(d);
    switchTab('molecules');
  });

  function updateMoleculeGallery(nodeData) {
    const grid = document.getElementById("mol-grid");
    const header = document.getElementById("mol-detail-header");
    header.innerText = `Molecules for ${nodeData.id} (${nodeData.run_id})`;
    
    if (!nodeData.molecules || nodeData.molecules.length === 0) {
      grid.innerHTML = `<div class="empty-state" style="grid-column: 1/-1; padding:40px">No molecules found for this node.</div>`;
      return;
    }

    grid.innerHTML = nodeData.molecules.map(m => `
      <div class="mol-card">
        <canvas class="mol-canvas" id="mol-${m.smiles.replace(/[^a-zA-Z0-9]/g,'')}"></canvas>
        <div class="mol-info">
          <div class="mol-be">${m.be != null ? m.be.toFixed(2) : 'N/A'} <span style="font-size:10px; font-weight:normal; color:var(--text-muted)">kcal/mol</span></div>
          <div class="mol-meta">
            <span>QED: ${m.qed != null ? m.qed.toFixed(2) : '?'}</span>
            <span>Steps: ${m.syn_steps || '?'}</span>
          </div>
          <div style="margin-top:8px">
            <span class="badge ${m.trivial ? 'badge-muted' : 'badge-success'}">
              ${m.trivial ? 'Trivial' : 'Valid Route'}
            </span>
          </div>
        </div>
      </div>
    `).join('');

    renderMolsWithRDKit(nodeData.molecules);
  }

  function renderMolsWithRDKit(mols) {
    if (window.RDKitModule) {
      mols.forEach(m => {
        const id = `mol-${m.smiles.replace(/[^a-zA-Z0-9]/g,'')}`;
        const canvas = document.getElementById(id);
        if (!canvas) return;
        try {
          const mol = RDKitModule.get_mol(m.smiles);
          mol.draw_to_canvas(canvas, canvas.width, canvas.height);
          mol.delete();
        } catch(e) { console.error("RDKit error", e); }
      });
    } else {
      const s = document.createElement("script");
      s.src = "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js";
      s.onload = () => {
        window.initRDKit().then(module => {
          window.RDKitModule = module;
          renderMolsWithRDKit(mols);
        });
      };
      document.head.appendChild(s);
    }
  }

  // Tab Navigation
  function switchTab(tabId) {
    document.querySelectorAll('.tab-btn').forEach(b => {
      b.classList.toggle('active', b.dataset.tab === tabId);
    });
    document.querySelectorAll('.tab-content').forEach(c => {
      c.classList.toggle('active', c.id === tabId);
    });
  }

  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => switchTab(btn.dataset.tab));
  });

  // Build History Table
  (function buildTable() {
    const hypoMap = {};
    GRAPH.filter(e => e.group === 'nodes' && e.data.id !== 'MOLCRAFT').forEach(el => {
      const id = el.data.id;
      if (!hypoMap[id]) hypoMap[id] = {};
      hypoMap[id][el.data.run_id] = el.data;
    });

    const hypoIds = Object.keys(hypoMap).sort();
    if (hypoIds.length === 0) return;

    let html = `<thead><tr><th>Hypothesis</th>`;
    SESSIONS.forEach(s => {
      html += `<th>${s.substring(5, 17)}</th>`;
    });
    html += `</tr></thead><tbody>`;

    hypoIds.forEach(hid => {
      html += `<tr><td style="font-weight:600">${hid}</td>`;
      let lastBe = null;
      SESSIONS.forEach(s => {
        const d = hypoMap[hid][s];
        if (d && d.best_be != null) {
          const be = d.best_be;
          let diffCls = "";
          if (lastBe !== null) {
            diffCls = be < lastBe ? "diff-up" : (be > lastBe ? "diff-down" : "");
          }
          const icon = d.success ? "✓" : "✗";
          html += `<td class="be-val ${diffCls}">${icon} ${be.toFixed(2)}</td>`;
          lastBe = be;
        } else {
          html += `<td style="color:var(--border)">-</td>`;
          lastBe = null;
        }
      });
      html += `</tr>`;
    });
    html += `</tbody>`;
    document.getElementById("session-table").innerHTML = html;
  })();

  // Legend
  const sessionColors = ["#38bdf8", "#fbbf24", "#34d399", "#a78bfa", "#f472b6", "#fb923c"];
  let legendHtml = SESSIONS.map((s, i) => `
    <div class="legend-item">
      <div class="dot" style="background:${sessionColors[i % sessionColors.length]}"></div>
      ${s.substring(5, 17)}
    </div>
  `).join('');
  legendHtml += `
    <div class="legend-item" style="margin-top:8px">
      <div class="dot" style="background:#22c55e"></div> Accepted
    </div>
    <div class="legend-item">
      <div class="dot" style="background:#475569"></div> Rejected
    </div>
  `;
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
    # Slate/Sky inspired palette
    session_colors = ["#38bdf8", "#fbbf24", "#34d399", "#a78bfa", "#f472b6", "#fb923c"]
    session_color_map = {s: session_colors[i % len(session_colors)] for i, s in enumerate(sessions)}

    for n in tree:
        hid = n["hypothesis_id"]
        success = n.get("success", False)
        run_id = n.get("run_id", "unknown")
        border = session_color_map.get(run_id, "#475569")
        # Accepted nodes are green, others are slate
        bg = "#22c55e" if success else "#475569"

        node = {
            "data": {
                "id": hid,
                "bg": bg,
                "border": border,
                "size": 32 if hid == "MOLCRAFT" else 28,
                "best_be": n.get("best_be"),
                "avg_be": n.get("avg_be"),
                "summary": n.get("summary", ""),
                "success": success,
                "run_id": run_id,
                "molecules": n.get("molecules", []),
            },
            "group": "nodes"
        }
        elements.append(node)

        parent = n.get("parent")
        if parent:
            edge = {
                "data": {
                    "id": f"{parent}_to_{hid}",
                    "source": parent,
                    "target": hid,
                    "ecolor": "#22c55e" if success else "#475569",
                    "estyle": "solid" if success else "dashed",
                },
                "group": "edges"
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
