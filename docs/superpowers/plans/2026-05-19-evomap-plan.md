# Molecular Evolution Map — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `tools/build_evomap.py` that reads `docs/iteration_log.jsonl` + `output/result.log` and produces a self-contained `output/evomap.html` with an interactive D3.js hypothesis evolution tree and RDKit JS molecule rendering.

**Architecture:** Single Python script (~180 lines) parses both JSONL sources, merges into a tree JSON, and embeds it into an inline HTML template. Output HTML loads D3.js v7 + RDKit JS from CDN. No pipeline modifications.

**Tech Stack:** Python 3 stdlib (json, pathlib, datetime), D3.js v7 CDN, RDKit JS minimallib CDN

---

## File Structure

| File | Action | Role |
|------|--------|------|
| `tools/build_evomap.py` | Create | Build script: parse, merge, embed, write |
| `tests/test_build_evomap.py` | Create | Unit tests for parse + merge logic |
| `output/evomap.html` | Generated | Self-contained deliverable |

No existing files modified.

---

## Data Reality Check

`iteration_log.jsonl` has 11 cumulative entries. `result.log` only has the most recent session's data. The tree structure (hypothesis lineage) comes from `iteration_log.jsonl`. Molecule details come from `result.log` for whatever rounds happen to match. Rounds without matching result.log data still get tree nodes — they just won't have molecule cards.

Result.log `start` events have a `round` field. Result.log `metrics` events have `min_binding_energy`, `avg_binding_energy`, `molecule_count`, `trivial_count`. Result.log `molecule` events have `mol_smiles`, `binding_energy`, `qed`, `syn_steps`, `trivial`, `route_quality`.

---

### Task 1: Data parsing functions

**Files:**
- Create: `tools/build_evomap.py` (partial — parse functions)
- Create: `tests/test_build_evomap.py` (partial — parse tests)

- [ ] **Step 1: Write failing tests for iteration_log parsing**

```python
import json
import tempfile
from pathlib import Path

def test_parse_iteration_log():
    tmp = Path(tempfile.mkstemp(suffix=".jsonl")[1])
    tmp.write_text(json.dumps({
        "round": 1, "hypothesis_id": "BASELINE", "success": True,
        "summary": "baseline run", "timestamp": "2026-05-17T15:39:04"
    }) + "\n" + json.dumps({
        "round": 1, "hypothesis_id": "H011", "success": True,
        "summary": "MMD diversity", "timestamp": "2026-05-18T05:36:05"
    }) + "\n")
    
    entries = parse_iteration_log(tmp)
    assert len(entries) == 2
    assert entries[0]["hypothesis_id"] == "BASELINE"
    assert entries[0]["best_be"] is None  # not in iteration_log
    assert entries[1]["summary"] == "MMD diversity"
    tmp.unlink()

def test_parse_result_log_start_and_metrics():
    tmp = Path(tempfile.mkstemp(suffix=".log")[1])
    tmp.write_text(json.dumps({
        "type": "start", "timestamp": "2026-05-18T21:19:51", "round": 1
    }) + "\n" + json.dumps({
        "type": "metrics", "timestamp": "2026-05-18T21:39:40",
        "molecule_count": 10, "trivial_count": 1,
        "avg_binding_energy": -8.06, "min_binding_energy": -9.17
    }) + "\n")
    
    runs = parse_result_log(tmp)
    assert len(runs) == 1
    assert runs[0]["round"] == 1
    assert runs[0]["metrics"]["min_binding_energy"] == -9.17
    assert runs[0]["metrics"]["trivial_ratio"] == 0.1
    tmp.unlink()

def test_parse_result_log_filters_non_events():
    tmp = Path(tempfile.mkstemp(suffix=".log")[1])
    tmp.write_text(json.dumps({
        "type": "stdout", "content": "[2026-05-18] Starting pipeline..."
    }) + "\n" + json.dumps({
        "type": "start", "timestamp": "2026-05-18T21:19:51", "round": 2
    }) + "\n")
    
    runs = parse_result_log(tmp)
    assert len(runs) == 1  # stdout event ignored
    tmp.unlink()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_build_evomap.py -v`
Expected: FAIL with "name 'parse_iteration_log' is not defined"

- [ ] **Step 3: Implement parse functions**

```python
"""build_evomap.py — 分子进化地图构建脚本。

读取 iteration_log.jsonl + result.log，生成自包含 D3.js HTML 进化树。
"""

import json
from pathlib import Path
from datetime import datetime
from typing import Any


def parse_iteration_log(path: Path) -> list[dict[str, Any]]:
    """Parse cumulative iteration log. One entry per Agent-reported round.

    Returns list of {round, hypothesis_id, success, summary, timestamp, best_be}.
    best_be is None here (not present in iteration_log); filled by merge step.
    """
    entries = []
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
                "best_be": None,    # filled by merge
                "avg_be": None,     # filled by merge
                "trivial_count": 0,
                "molecule_count": 0,
                "molecules": [],    # filled by merge
            })
    return entries


def parse_result_log(path: Path) -> list[dict[str, Any]]:
    """Parse structured result.log. Groups events by run (delimited by 'start' events).

    Returns list of runs, each: {round, start_ts, metrics: {...}, molecules: [...]}.
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
                # New run delimiter
                if current_run is not None and current_run.get("molecules"):
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

    # Don't forget the last run
    if current_run is not None and current_run.get("molecules"):
        runs.append(current_run)

    return runs
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_build_evomap.py::test_parse_iteration_log tests/test_build_evomap.py::test_parse_result_log_start_and_metrics tests/test_build_evomap.py::test_parse_result_log_filters_non_events -v`
Expected: 3 PASS

- [ ] **Step 5: Commit**

```bash
git add tools/build_evomap.py tests/test_build_evomap.py
git commit -m "feat: add data parsing functions for evomap build script"
```

---

### Task 2: Tree construction and merge logic

**Files:**
- Modify: `tools/build_evomap.py` (add merge + tree functions)
- Modify: `tests/test_build_evomap.py` (add merge tests)

- [ ] **Step 1: Write failing test for merge**

```python
def test_merge_builds_tree():
    iter_entries = [
        {"round": 1, "hypothesis_id": "BASELINE", "success": True,
         "summary": "base", "timestamp": "2026-05-17T10:00:00",
         "best_be": None, "avg_be": None, "trivial_count": 0, "molecule_count": 0, "molecules": []},
        {"round": 1, "hypothesis_id": "H001", "success": True,
         "summary": "first hypothesis", "timestamp": "2026-05-17T11:00:00",
         "best_be": None, "avg_be": None, "trivial_count": 0, "molecule_count": 0, "molecules": []},
        {"round": 2, "hypothesis_id": "H002", "success": False,
         "summary": "failed one", "timestamp": "2026-05-17T12:00:00",
         "best_be": None, "avg_be": None, "trivial_count": 0, "molecule_count": 0, "molecules": []},
    ]
    runs = [
        {"round": 1, "start_ts": "2026-05-17T10:00:00",
         "metrics": {"min_binding_energy": -8.56, "avg_binding_energy": -8.06,
                     "molecule_count": 10, "trivial_count": 0, "trivial_ratio": 0.0},
         "molecules": [
             {"smiles": "c1ccccc1", "be": -8.56, "qed": 0.8, "trivial": False, "syn_steps": 1}
         ]},
    ]
    
    tree = build_tree(iter_entries, runs)
    assert len(tree) == 3
    # First entry should have metrics merged from matching run
    assert tree[0]["best_be"] == -8.56
    assert tree[0]["avg_be"] == -8.06
    assert len(tree[0]["molecules"]) == 1
    # Second entry has no matching run data
    assert tree[1]["best_be"] is None
    # All entries should have parent assigned
    assert tree[0]["parent"] is None  # root
    assert tree[1]["parent"] == "BASELINE"
    assert tree[2]["parent"] == "H001"

def test_merge_no_result_log():
    """Tree should still be built even when result.log is empty."""
    iter_entries = [
        {"round": 1, "hypothesis_id": "BASELINE", "success": True,
         "summary": "base", "timestamp": "2026-05-17T10:00:00",
         "best_be": None, "avg_be": None, "trivial_count": 0, "molecule_count": 0, "molecules": []},
    ]
    tree = build_tree(iter_entries, [])
    assert len(tree) == 1
    assert tree[0]["best_be"] is None  # gracefully missing
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest tests/test_build_evomap.py::test_merge_builds_tree tests/test_build_evomap.py::test_merge_no_result_log -v`
Expected: FAIL with "name 'build_tree' is not defined"

- [ ] **Step 3: Implement merge + tree functions**

```python
def _extract_best_be(summary: str) -> float | None:
    """Try to extract best BE from summary text if available."""
    import re
    # Match patterns like "最佳结合能 -8.56" or "best BE -8.56" or "-8.56 kcal/mol"
    m = re.search(r'([-]?\d+\.\d+)\s*kcal/mol', summary)
    if m:
        return float(m.group(1))
    return None


def build_tree(
    iter_entries: list[dict[str, Any]],
    runs: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Merge iteration log entries with run metrics, assign parentage.

    Matching strategy: align by chronological order — the N-th run matches
    the N-th iteration entry with the same round number. If there's only
    one run in result.log, its metrics go to the LAST matching round entry.

    Parentage: entries are in chronological order. Each entry's parent is
    the most recent ACCEPTED entry before it. Root has parent=None.
    """
    # Build run lookup by round → list of runs (preserving order)
    run_by_round: dict[int, list[dict[str, Any]]] = {}
    for run in runs:
        rn = run["round"]
        run_by_round.setdefault(rn, []).append(run)

    # Merge: for each iter entry, try to find matching run
    run_consumed: dict[int, int] = {}  # round → index into run_by_round list

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
            # Only keep top-5 molecules by BE for display
            mols = sorted(run.get("molecules", []),
                         key=lambda x: x.get("be") or 999)
            entry["molecules"] = mols[:5]
        else:
            # No matching run — try to extract BE from summary text
            entry["best_be"] = _extract_best_be(entry["summary"])

    # Assign parentage: parent = most recent prior ACCEPTED hypothesis
    last_accepted: dict[int, str] = {}  # round → hypothesis_id
    for entry in iter_entries:
        rn = entry["round"]
        # Find parent: accepted entry from an earlier round
        parent = None
        for pr in sorted(last_accepted.keys(), reverse=True):
            if pr < rn or (pr == rn and last_accepted[pr] != entry["hypothesis_id"]):
                parent = last_accepted[pr]
                break
        entry["parent"] = parent
        if entry["success"]:
            last_accepted[rn] = entry["hypothesis_id"]

    return iter_entries
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_build_evomap.py::test_merge_builds_tree tests/test_build_evomap.py::test_merge_no_result_log -v`
Expected: 2 PASS

- [ ] **Step 5: Commit**

```bash
git add tools/build_evomap.py tests/test_build_evomap.py
git commit -m "feat: add tree construction and merge logic"
```

---

### Task 3: HTML generation with D3.js tree + RDKit JS

**Files:**
- Modify: `tools/build_evomap.py` (add HTML template + embed function)

- [ ] **Step 1: Write failing test for HTML generation**

```python
def test_generate_html_writes_valid_file():
    import tempfile
    tree = [{
        "id": "BASELINE", "round": 1, "success": True,
        "hypothesis_id": "BASELINE", "parent": None,
        "best_be": -8.56, "avg_be": -8.06, "trivial_count": 0,
        "molecule_count": 10, "summary": "baseline",
        "molecules": [
            {"smiles": "c1ccccc1", "be": -8.56, "qed": 0.72, "trivial": False, "syn_steps": 1}
        ]
    }]
    
    out = Path(tempfile.mkstemp(suffix=".html")[1])
    generate_html(tree, out)
    
    html = out.read_text()
    assert "<!DOCTYPE html>" in html
    assert "BASELINE" in html
    assert "d3.js" in html.lower() or "d3@" in html
    assert "rdkit" in html.lower() or "RDKit" in html
    assert "c1ccccc1" in html  # molecule data embedded
    out.unlink()

def test_generate_html_empty_tree():
    import tempfile
    out = Path(tempfile.mkstemp(suffix=".html")[1])
    generate_html([], out)
    html = out.read_text()
    assert "<!DOCTYPE html>" in html
    assert "暂无数据" in html or "No data" in html  # graceful empty state
    out.unlink()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest tests/test_build_evomap.py::test_generate_html_writes_valid_file tests/test_build_evomap.py::test_generate_html_empty_tree -v`
Expected: FAIL with "name 'generate_html' is not defined"

- [ ] **Step 3: Implement HTML template and generate function**

```python
HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Molecular Evolution Map — MolCraft Agent</title>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; background: #1a1a2e; color: #eee; overflow: hidden; height: 100vh; }
#tree-container { width: 100%; height: 100%; position: relative; }

/* Node styles */
.node circle { stroke-width: 2px; }
.node .success circle { fill: #2ecc71; stroke: #27ae60; }
.node .failure circle { fill: #95a5a6; stroke: #7f8c8d; }
.node text { font-size: 12px; fill: #ecf0f1; font-family: monospace; }
.node .be-label { font-size: 10px; fill: #bdc3c7; }

/* Edge styles */
.link { fill: none; stroke-width: 2px; }
.link.success { stroke: #27ae60; }
.link.failure { stroke: #7f8c8d; stroke-dasharray: 6,4; }

/* Tooltip */
.tooltip {
  position: absolute; padding: 10px 14px; background: rgba(44,62,80,0.95);
  border: 1px solid #34495e; border-radius: 6px; font-size: 13px;
  pointer-events: none; max-width: 360px; line-height: 1.5;
  box-shadow: 0 4px 12px rgba(0,0,0,0.4);
}
.tooltip .hyp-id { font-weight: bold; color: #3498db; font-size: 15px; }
.tooltip .be-big { color: #e74c3c; font-size: 16px; font-weight: bold; }

/* Detail panel */
#detail-panel {
  position: absolute; bottom: 0; left: 0; right: 0; height: 200px;
  background: rgba(30,39,56,0.97); border-top: 1px solid #34495e;
  overflow-x: auto; overflow-y: hidden; white-space: nowrap;
  padding: 12px 16px; display: none;
}
#detail-panel .mol-card {
  display: inline-block; width: 180px; margin-right: 12px;
  background: #243447; border-radius: 6px; padding: 8px;
  text-align: center; vertical-align: top;
}
.mol-card .be { color: #e74c3c; font-weight: bold; }
.mol-card .qed { color: #bdc3c7; font-size: 11px; }
.mol-card .route-badge {
  display: inline-block; padding: 2px 6px; border-radius: 3px;
  font-size: 10px; margin-top: 4px;
}
.mol-card .route-badge.trivial { background: #7f8c8d; color: #fff; }
.mol-card .route-badge.non-trivial { background: #27ae60; color: #fff; }

.empty-state { text-align: center; padding: 40px; color: #7f8c8d; font-size: 16px; }
</style>
</head>
<body>
<div id="tree-container">
  <svg width="100%" height="100%"></svg>
  <div class="tooltip" style="display:none"></div>
</div>
<div id="detail-panel"><div style="padding:8px;color:#bdc3c7;font-size:12px;">点击节点查看分子详情</div></div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
document.addEventListener("DOMContentLoaded", () => {
  const DATA = __DATA_PLACEHOLDER__;

  if (!DATA || DATA.length === 0) {
    document.getElementById("tree-container").innerHTML =
      '<div class="empty-state">暂无进化数据<br><small>运行 main.py 完成至少一轮迭代后执行 build_evomap.py</small></div>';
    return;
  }

  const svg = d3.select("#tree-container svg");
  const tooltip = d3.select(".tooltip");
  const detailPanel = d3.select("#detail-panel");
  const margin = {top: 40, right: 120, bottom: 40, left: 160};

  // Build D3 hierarchy
  const rootData = DATA.find(d => d.parent === null) || DATA[0];
  const stratify = d3.stratify()
    .id(d => d.hypothesis_id)
    .parentId(d => d.parent);
  const root = stratify(DATA);

  // Assign node dimensions
  root.each(d => { d._children = d.children; });

  const width = window.innerWidth;
  const height = window.innerHeight - 200;
  const treeLayout = d3.tree().size([height - margin.top - margin.bottom, width - margin.left - margin.right]);
  treeLayout(root);

  // Draw edges
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

  g.selectAll(".link")
    .data(root.links())
    .join("path")
    .attr("class", d => `link ${d.target.data.success ? "success" : "failure"}`)
    .attr("d", d3.linkHorizontal()
      .x(d => d.y)
      .y(d => d.x))
    .append("title")
    .text(d => `best BE: ${d.target.data.best_be || "N/A"}`);

  // Draw nodes
  const node = g.selectAll(".node")
    .data(root.descendants())
    .join("g")
    .attr("class", d => `node ${d.data.success ? "success" : "failure"}`)
    .attr("transform", d => `translate(${d.y},${d.x})`)
    .on("mouseover", (event, d) => {
      const data = d.data;
      tooltip.style("display", "block")
        .html(`<div class="hyp-id">${data.hypothesis_id}</div>
               <div>最佳 BE: <span class="be-big">${data.best_be != null ? data.best_be.toFixed(3) + " kcal/mol" : "N/A"}</span></div>
               <div>平均 BE: ${data.avg_be != null ? data.avg_be.toFixed(3) + " kcal/mol" : "N/A"}</div>
               <div>Molecules: ${data.molecule_count || "N/A"} | Trivial: ${data.trivial_count}</div>
               <div style="margin-top:6px;color:#bdc3c7;font-size:11px;max-height:60px;overflow:hidden;">${data.summary.substring(0, 200)}</div>`);
    })
    .on("mousemove", (event) => {
      tooltip.style("left", (event.pageX + 12) + "px")
        .style("top", (event.pageY - 28) + "px");
    })
    .on("mouseout", () => tooltip.style("display", "none"))
    .on("click", (event, d) => {
      showMolecules(d.data);
    });

  node.append("circle")
    .attr("r", d => d.data.hypothesis_id === "BASELINE" ? 8 : 6);

  node.append("text")
    .attr("dy", -12)
    .attr("text-anchor", "middle")
    .text(d => d.data.hypothesis_id);

  node.filter(d => d.data.best_be != null).append("text")
    .attr("dy", 20)
    .attr("text-anchor", "middle")
    .attr("class", "be-label")
    .text(d => d.data.best_be.toFixed(1));

  function showMolecules(data) {
    const mols = data.molecules || [];
    if (mols.length === 0) {
      detailPanel.style("display", "block")
        .html(`<div style="padding:12px;color:#7f8c8d;">${data.hypothesis_id}: 本轮没有保存分子数据</div>`);
      return;
    }

    let cards = `<div style="padding:4px 0 4px 12px;color:#bdc3c7;font-size:12px;">${data.hypothesis_id} — ${data.summary.substring(0, 120)}</div>`;
    mols.forEach(m => {
      const isTriv = m.trivial;
      cards += `<div class="mol-card">
        <div id="mol-${m.smiles.replace(/[^a-zA-Z0-9]/g,'')}" style="width:160px;height:100px;margin:0 auto;"></div>
        <div class="be">BE: ${m.be != null ? m.be.toFixed(2) : "N/A"}</div>
        <div class="qed">QED: ${m.qed != null ? m.qed.toFixed(2) : "N/A"} | Steps: ${m.syn_steps || "?"}</div>
        <span class="route-badge ${isTriv ? 'trivial' : 'non-trivial'}">${isTriv ? 'trivial' : 'valid route'}</span>
      </div>`;
    });
    detailPanel.style("display", "block").html(cards);

    // Render molecules via RDKit JS
    if (typeof initRDKit === "undefined") {
      // RDKit not loaded yet — inject script
      const script = document.createElement("script");
      script.src = "https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js";
      script.onload = () => initRDKit().then(() => renderMolCards(mols));
      document.head.appendChild(script);
    } else {
      initRDKit().then(() => renderMolCards(mols));
    }
  }

  function renderMolCards(mols) {
    mols.forEach(m => {
      const divId = "mol-" + m.smiles.replace(/[^a-zA-Z0-9]/g, '');
      const el = document.getElementById(divId);
      if (!el) return;
      try {
        const mol = RDKitModule.get_mol(m.smiles);
        if (!mol) return;
        mol.draw_to_canvas(el, 160, 100);
        mol.delete();
      } catch(e) {
        el.innerHTML = '<span style="color:#e74c3c;font-size:10px;">render error</span>';
      }
    });
  }

  // Handle resize
  window.addEventListener("resize", () => {
    const w = window.innerWidth, h = window.innerHeight;
    svg.attr("width", w).attr("height", h);
  });
});
</script>
</body>
</html>
'''


def generate_html(tree: list[dict[str, Any]], output_path: Path) -> None:
    """Embed tree data into HTML template and write to output_path.

    Uses parse_constant to reject NaN/Infinity (same pattern as main.py).
    """
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_build_evomap.py::test_generate_html_writes_valid_file tests/test_build_evomap.py::test_generate_html_empty_tree -v`
Expected: 2 PASS

- [ ] **Step 5: Commit**

```bash
git add tools/build_evomap.py tests/test_build_evomap.py
git commit -m "feat: add D3.js HTML generation with molecule rendering"
```

---

### Task 4: CLI entry point and integration test

**Files:**
- Modify: `tools/build_evomap.py` (add main function)
- Modify: `tests/test_build_evomap.py` (add integration test)

- [ ] **Step 1: Write failing integration test**

```python
def test_end_to_end_with_real_files():
    """Test against actual project data files."""
    import subprocess, sys
    project_root = Path(__file__).parent.parent
    iter_log = project_root / "docs" / "iteration_log.jsonl"
    result_log = project_root / "output" / "result.log"
    
    if not iter_log.exists():
        pytest.skip("iteration_log.jsonl not found")
    
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "evomap.html"
        result = subprocess.run(
            [sys.executable, str(project_root / "tools" / "build_evomap.py"),
             "--iteration-log", str(iter_log),
             "--result-log", str(result_log) if result_log.exists() else "",
             "--output", str(out)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"Build failed: {result.stderr}"
        assert out.exists()
        html = out.read_text()
        assert "BASELINE" in html
        assert "H011" in html
        assert "H015" in html  # latest hypothesis
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest tests/test_build_evomap.py::test_end_to_end_with_real_files -v`
Expected: FAIL (no main function / argparse)

- [ ] **Step 3: Implement main function**

```python
def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Build Molecular Evolution Map HTML"
    )
    parser.add_argument(
        "--iteration-log",
        type=Path,
        default=Path("docs/iteration_log.jsonl"),
        help="Path to iteration_log.jsonl (default: docs/iteration_log.jsonl)",
    )
    parser.add_argument(
        "--result-log",
        type=Path,
        default=Path("output/result.log"),
        help="Path to result.log (default: output/result.log)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/evomap.html"),
        help="Output HTML path (default: output/evomap.html)",
    )
    args = parser.parse_args()

    # Resolve relative paths from project root
    project_root = Path(__file__).resolve().parent.parent
    iter_log = args.iteration_log if args.iteration_log.is_absolute() else project_root / args.iteration_log
    result_log = args.result_log if args.result_log.is_absolute() else project_root / args.result_log
    output_path = args.output if args.output.is_absolute() else project_root / args.output

    if not iter_log.exists():
        print(f"✗ 未找到 iteration_log: {iter_log}")
        print("  请先运行 main.py 完成至少一轮迭代。")
        raise SystemExit(1)

    print(f"读取迭代日志: {iter_log}")
    iter_entries = parse_iteration_log(iter_log)
    print(f"  共 {len(iter_entries)} 条迭代记录")

    runs = []
    if result_log.exists():
        print(f"读取运行日志: {result_log}")
        runs = parse_result_log(result_log)
        print(f"  共 {len(runs)} 次 pipeline 运行")
    else:
        print(f"未找到 result.log，树节点将不包含分子详情")

    print("构建进化树...")
    tree = build_tree(iter_entries, runs)

    print(f"生成 HTML: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    generate_html(tree, output_path)

    print(f"\n✓ 完成! 在浏览器中打开: {output_path}")
    print(f"  节点数: {len(tree)}")
    n_success = sum(1 for n in tree if n["success"])
    n_fail = sum(1 for n in tree if not n["success"])
    print(f"  成功: {n_success} | 失败: {n_fail}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run integration test**

Run: `python3 -m pytest tests/test_build_evomap.py::test_end_to_end_with_real_files -v`
Expected: PASS

- [ ] **Step 5: Run against real data and verify visually**

Run: `python3 tools/build_evomap.py`
Expected: `output/evomap.html` created, open in browser to verify tree renders correctly

- [ ] **Step 6: Commit**

```bash
git add tools/build_evomap.py tests/test_build_evomap.py
git commit -m "feat: add CLI entry point for evomap build script"
```
