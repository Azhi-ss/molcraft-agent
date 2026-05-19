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
