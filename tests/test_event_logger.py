"""Tests for EventLogger — output format, atomic commit, edge cases."""

import json
import os
import tempfile
from pathlib import Path

import pytest

from src.event_logger import EventLogger
from src.event_schema import MoleculeEvent, MetricsEvent


@pytest.fixture
def log_path():
    tmpdir = tempfile.mkdtemp()
    yield Path(tmpdir) / "result.log"
    # cleanup
    for p in Path(tmpdir).iterdir():
        p.unlink(missing_ok=True)
    os.rmdir(tmpdir)


@pytest.fixture
def logger(log_path):
    """Provide a logger that auto-cleans up."""
    lgr = EventLogger(log_path)
    yield lgr
    lgr.close()


def _read_log_content(log_path: Path) -> list[dict]:
    with open(log_path) as f:
        return [json.loads(line) for line in f if line.strip()]


class TestEventLoggerOutput:
    def test_log_start_writes_valid_jsonl(self, logger, log_path):
        logger.log_start(round=1, max_minutes=120, max_iterations=3, program_file="program.md")
        logger.commit()
        events = _read_log_content(log_path)
        assert len(events) == 1
        assert events[0]["type"] == "start"
        assert events[0]["round"] == 1
        assert "timestamp" in events[0]

    def test_log_end_writes_status(self, logger, log_path):
        logger.log_end(status="success")
        logger.commit()
        events = _read_log_content(log_path)
        assert events[0]["type"] == "end"
        assert events[0]["status"] == "success"

    def test_log_stage_begun(self, logger, log_path):
        logger.log_stage(name="诊断", status="begun")
        logger.commit()
        events = _read_log_content(log_path)
        assert events[0]["type"] == "stage"
        assert events[0]["name"] == "诊断"
        assert events[0]["status"] == "begun"
        assert "duration_seconds" not in events[0]

    def test_log_stage_completed(self, logger, log_path):
        logger.log_stage(name="诊断", status="completed", duration_seconds=318.0)
        logger.commit()
        events = _read_log_content(log_path)
        assert events[0]["duration_seconds"] == 318.0

    def test_log_molecule_writes_all_fields(self, logger, log_path):
        ev = MoleculeEvent(
            mol_smiles="Cc1ccccc1",
            binding_energy=-8.5,
            composite_score=0.97,
            syn_steps=2,
            sa_score=3.2,
            qed=0.68,
            rings=2,
            docking_std=0.05,
            route_quality=0.85,
        )
        logger.log_molecule(ev)
        logger.commit()
        events = _read_log_content(log_path)
        assert events[0]["type"] == "molecule"
        assert events[0]["mol_smiles"] == "Cc1ccccc1"
        assert events[0]["binding_energy"] == -8.5
        assert events[0]["composite_score"] == 0.97
        assert events[0]["syn_steps"] == 2
        assert events[0]["sa_score"] == 3.2
        assert events[0]["qed"] == 0.68
        assert events[0]["rings"] == 2
        assert events[0]["docking_std"] == 0.05
        assert events[0]["route_quality"] == 0.85

    def test_log_metrics_unified_shape(self, logger, log_path):
        ev = MetricsEvent(
            molecule_count=10,
            non_trivial_count=8,
            trivial_count=2,
            avg_binding_energy=-7.5,
            min_binding_energy=-9.0,
            avg_syn_steps=1.5,
            docking_success_rate=0.9,
        )
        logger.log_metrics(ev)
        logger.commit()
        events = _read_log_content(log_path)
        assert events[0]["type"] == "metrics"
        assert events[0]["molecule_count"] == 10
        assert events[0]["non_trivial_count"] == 8
        assert events[0]["trivial_count"] == 2
        assert events[0]["docking_success_rate"] == 0.9

    def test_log_hypothesis_validation(self, logger, log_path):
        logger.log_hypothesis_validation(
            hypothesis_id="H001", success=True, conclusion="ACCEPTED",
            changes={"threshold": {"from": 8, "to": 6}},
        )
        logger.commit()
        events = _read_log_content(log_path)
        assert events[0]["type"] == "hypothesis_validation"
        assert events[0]["hypothesis_id"] == "H001"
        assert events[0]["success"] is True
        assert events[0]["changes"]["threshold"]["from"] == 8


class TestEventLoggerLifecycle:
    def test_log_writes_to_temp_before_commit(self, logger, log_path):
        logger.log_start(round=1, max_minutes=120, max_iterations=3, program_file="p.md")
        assert not log_path.exists()
        assert logger.tmp_path.exists()

    def test_commit_creates_final_file(self, logger, log_path):
        logger.log_start(round=1, max_minutes=120, max_iterations=3, program_file="p.md")
        logger.commit()
        assert log_path.exists()
        assert not logger.tmp_path.exists()

    def test_close_without_commit_removes_temp(self, logger, log_path):
        logger.log_start(round=1, max_minutes=120, max_iterations=3, program_file="p.md")
        tmp = logger.tmp_path
        logger.close()
        assert not tmp.exists()
        assert not log_path.exists()

    def test_stdout_writes_event(self, logger, log_path):
        logger.write("Hello world\n")
        logger.commit()
        events = _read_log_content(log_path)
        assert len(events) == 1
        assert events[0]["type"] == "stdout"
        assert "Hello world" in events[0]["content"]


class TestTruncation:
    """Test that verbose agent output is truncated in log but not terminal."""

    def test_file_read_dump_truncated(self, logger, log_path):
        """ReadFile output: keep only system line + first 3 content lines."""
        msg = (
            "<system>300 lines read from file starting from line 1. End of file reached.</system>\n"
            + "\n".join(f"     {i}\tcontent line {i}" for i in range(1, 301))
        )
        logger.write(msg)
        logger.commit()
        events = _read_log_content(log_path)
        content = events[0]["content"]
        assert "lines read from file" in content
        assert "content line 1" in content
        assert "content line 3" in content
        assert "content line 4" not in content  # truncated after line 3
        assert "truncated" in content.lower()
        assert len(content) < 1000  # was thousands of chars

    def test_normal_output_passes_through(self, logger, log_path):
        """Short agent messages are NOT truncated."""
        msg = "Pipeline completed: 10 molecules, best BE=-10.1 kcal/mol\n"
        logger.write(msg)
        logger.commit()
        events = _read_log_content(log_path)
        assert "Pipeline completed" in events[0]["content"]
        assert "truncated" not in events[0]["content"].lower()

    def test_code_file_read_truncated(self, logger, log_path):
        """ReadFile of long code files is also truncated."""
        lines = ["<system>200 lines read from file. End of file reached.</system>"]
        lines.extend(f"     {i}\tdef long_function_{i}(): pass  # line {i}" for i in range(1, 201))
        msg = "\n".join(lines)
        logger.write(msg)
        logger.commit()
        events = _read_log_content(log_path)
        content = events[0]["content"]
        assert "lines read from file" in content
        assert "long_function_1" in content
        assert "long_function_4" not in content  # truncated

    def test_generic_long_output_truncated(self, logger, log_path):
        """Any stdout line > 2000 chars that isn't ReadFile gets hard-truncated."""
        msg = "x" * 3000
        logger.write(msg)
        logger.commit()
        events = _read_log_content(log_path)
        assert len(events[0]["content"]) <= 2100  # 2000 + "...[truncated]"


class TestEdgeCases:
    def test_nan_rejected(self, logger, log_path):
        ev = MoleculeEvent(
            mol_smiles="Cc1ccccc1",
            binding_energy=float("nan"),
            composite_score=0.97,
            syn_steps=2,
            sa_score=3.2,
            qed=0.68,
            rings=2,
        )
        with pytest.raises(ValueError):
            logger.log_molecule(ev)

    def test_empty_string_write_no_event(self, logger, log_path):
        logger.write("")
        logger.write("\n")
        logger.write("   ")
        logger.commit()
        with open(log_path) as f:
            content = f.read().strip()
        assert content == ""  # no stdout events for blank content

    def test_multiple_events_sequential(self, logger, log_path):
        logger.log_start(round=1, max_minutes=120, max_iterations=3, program_file="p.md")
        logger.log_end(status="success")
        logger.commit()
        events = _read_log_content(log_path)
        assert len(events) == 2
        assert events[0]["type"] == "start"
        assert events[1]["type"] == "end"
