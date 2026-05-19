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
) -> list[dict[str, Any]]:
    """Merge iteration log entries with run metrics and assign parentage.

    Matching: runs aligned to iter entries by round number + chronological order.
    Parentage: each entry's parent is the most recent ACCEPTED hypothesis before it.
    Root entries have parent=None.
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
