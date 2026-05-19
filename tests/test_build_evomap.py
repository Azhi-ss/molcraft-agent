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


# --- Task 2: Tree construction tests ---


def test_merge_builds_tree():
    from tools.build_evomap import build_tree

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
    assert tree[0]["best_be"] == -8.56
    assert tree[0]["avg_be"] == -8.06
    assert len(tree[0]["molecules"]) == 1
    assert tree[1]["best_be"] is None  # no matching run
    assert tree[0]["parent"] is None   # root
    assert tree[1]["parent"] == "BASELINE"
    assert tree[2]["parent"] == "H001"


def test_merge_no_result_log():
    from tools.build_evomap import build_tree

    iter_entries = [
        {"round": 1, "hypothesis_id": "BASELINE", "success": True,
         "summary": "base", "timestamp": "2026-05-17T10:00:00",
         "best_be": None, "avg_be": None, "trivial_count": 0, "molecule_count": 0, "molecules": []},
    ]
    tree = build_tree(iter_entries, [])
    assert len(tree) == 1
    assert tree[0]["best_be"] is None


def test_extract_best_be_from_summary():
    from tools.build_evomap import _extract_best_be
    assert _extract_best_be("最佳结合能从 -8.335 提升至 -9.941") == -9.941
    assert _extract_best_be("best BE -8.56 kcal/mol") == -8.56
    assert _extract_best_be("no energy here") is None
    assert _extract_best_be("") is None


# --- Task 3: HTML generation tests ---


def test_generate_html_writes_valid_file():
    from tools.build_evomap import generate_html

    tree = [{
        "hypothesis_id": "BASELINE", "round": 1, "success": True,
        "parent": None, "best_be": -8.56, "avg_be": -8.06,
        "trivial_count": 0, "molecule_count": 10,
        "summary": "baseline",
        "molecules": [
            {"smiles": "c1ccccc1", "be": -8.56, "qed": 0.72,
             "trivial": False, "syn_steps": 1}
        ]
    }]

    out = Path(tempfile.mkstemp(suffix=".html")[1])
    generate_html(tree, out)

    html = out.read_text()
    assert "<!DOCTYPE html>" in html
    assert "BASELINE" in html
    assert "cytoscape.min.js" in html
    assert "RDKit" in html or "rdkit" in html
    assert "c1ccccc1" in html
    out.unlink()


def test_generate_html_empty_tree():
    from tools.build_evomap import generate_html

    out = Path(tempfile.mkstemp(suffix=".html")[1])
    generate_html([], out)
    html = out.read_text()
    assert "<!DOCTYPE html>" in html
    out.unlink()


# --- Molecules JSONL test ---


def test_parse_molecules_jsonl():
    from tools.build_evomap import parse_molecules_jsonl

    tmp = Path(tempfile.mkstemp(suffix=".jsonl")[1])
    tmp.write_text(json.dumps({
        "timestamp": "2026-05-18T21:19:51",
        "molecules": [
            {"smiles": "c1ccccc1", "be": -8.56, "qed": 0.8, "trivial": False},
            {"smiles": "c1ccccc1O", "be": -8.10, "qed": 0.7, "trivial": True},
        ]
    }) + "\n")

    runs = parse_molecules_jsonl(tmp)
    assert len(runs) == 1
    assert len(runs[0]["molecules"]) == 2
    assert runs[0]["molecules"][0]["be"] == -8.56
    tmp.unlink()


def test_build_tree_with_molecules():
    from tools.build_evomap import build_tree

    iter_entries = [
        {"round": 1, "hypothesis_id": "BASELINE", "success": True,
         "summary": "base", "timestamp": "2026-05-17T10:00:00",
         "best_be": None, "avg_be": None, "trivial_count": 0, "molecule_count": 0, "molecules": []},
    ]
    mol_runs = [
        {"timestamp": "2026-05-17T11:00:00",
         "molecules": [
             {"smiles": "c1ncccc1", "be": -9.0, "qed": 0.9, "trivial": False}
         ]},
    ]

    tree = build_tree(iter_entries, [], mol_runs)
    assert len(tree) == 1
    assert len(tree[0]["molecules"]) == 1
    assert tree[0]["molecules"][0]["smiles"] == "c1ncccc1"


# --- Task 4: Integration test ---


def test_end_to_end_with_real_files():
    """Build evomap from actual project data. Skips if sources missing."""
    from tools.build_evomap import parse_iteration_log, parse_result_log, build_tree, generate_html

    project_root = Path(__file__).parent.parent
    iter_log = project_root / "docs" / "iteration_log.jsonl"
    result_log = project_root / "output" / "result.log"

    if not iter_log.exists():
        pytest.skip("iteration_log.jsonl not found")

    entries = parse_iteration_log(iter_log)
    assert len(entries) >= 1

    runs = parse_result_log(result_log) if result_log.exists() else []
    tree = build_tree(entries, runs)
    assert len(tree) >= 1
    assert tree[0]["hypothesis_id"] is not None

    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "evomap.html"
        generate_html(tree, out)
        assert out.exists()
        html = out.read_text()
        assert "H015" in html or "BASELINE" in html
