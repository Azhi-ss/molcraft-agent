import json
import tempfile
import subprocess
import sys
from pathlib import Path

import pytest

# --- Task 1: Data parsing tests ---


def test_parse_iteration_log():
    from tools.build_evomap import parse_iteration_log

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
    assert entries[0]["best_be"] is None
    assert entries[1]["summary"] == "MMD diversity"
    tmp.unlink()


def test_parse_result_log_start_and_metrics():
    from tools.build_evomap import parse_result_log

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
    from tools.build_evomap import parse_result_log

    tmp = Path(tempfile.mkstemp(suffix=".log")[1])
    tmp.write_text(json.dumps({
        "type": "stdout", "content": "[2026-05-18] Starting pipeline..."
    }) + "\n" + json.dumps({
        "type": "start", "timestamp": "2026-05-18T21:19:51", "round": 2
    }) + "\n")

    runs = parse_result_log(tmp)
    assert len(runs) == 1
    tmp.unlink()
