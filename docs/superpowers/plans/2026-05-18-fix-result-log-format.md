# result.log JSON Format Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract EventLogger module, unify metrics fields, add molecule/stage events, TDD-driven.

**Architecture:** New `src/event_schema.py` defines dataclasses with validation. New `src/event_logger.py` provides EventLogger (replaces StructuredLogger). `main.py` and `pipeline.py` consume it; `tools.py` adds Agent-triggered stage tools.

**Tech Stack:** Python 3.10+, dataclasses, json, pytest

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `src/event_schema.py` | Dataclass definitions per event type, `to_dict()`, field validation |
| Create | `src/event_logger.py` | EventLogger class: JSONL write, atomic commit, terminal passthrough |
| Create | `tests/test_event_schema.py` | Schema validation tests |
| Create | `tests/test_event_logger.py` | EventLogger unit tests |
| Modify | `main.py` | Replace StructuredLogger → EventLogger, remove pipeline log merge |
| Modify | `tools/pipeline.py` | Remove `write_result_log()`, receive EventLogger, call `log_molecule()` per candidate |
| Modify | `tools/tools.py` | Add `begin_stage(name)` / `end_stage()` tools |
| Modify | `program.md` | Instruct Agent to call begin_stage/end_stage at phase transitions |

---

### Task 1: event dataclass definitions

**Files:**
- Create: `src/event_schema.py`
- Test: (tested in Task 2)

- [ ] **Step 1: Write `src/event_schema.py`**

```python
"""Event type dataclasses for structured JSONL logging.

Each event type is a dataclass with a to_dict() method for JSON serialization.
Validation occurs at construction time.
"""

from dataclasses import dataclass, field, asdict
from typing import Any


def _check_required(value: Any, field_name: str) -> None:
    if value is None:
        raise ValueError(f"{field_name} is required")


@dataclass
class StartEvent:
    round: int
    max_minutes: int
    max_iterations: int
    program_file: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class EndEvent:
    status: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class StageEvent:
    name: str
    status: str  # "begun" | "completed" | "failed"
    duration_seconds: float | None = None

    def __post_init__(self) -> None:
        valid_statuses = ("begun", "completed", "failed")
        if self.status not in valid_statuses:
            raise ValueError(f"status must be one of {valid_statuses}, got {self.status!r}")

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        # Remove None values for cleaner output
        if d.get("duration_seconds") is None:
            del d["duration_seconds"]
        return d


@dataclass
class MoleculeEvent:
    mol_smiles: str
    binding_energy: float
    composite_score: float
    syn_steps: int
    sa_score: float | None
    qed: float
    rings: int
    trivial: bool = False
    route_quality: float | None = None
    docking_std: float | None = None

    def __post_init__(self) -> None:
        _check_required(self.mol_smiles, "mol_smiles")
        _check_required(self.binding_energy, "binding_energy")
        _check_required(self.composite_score, "composite_score")
        _check_required(self.syn_steps, "syn_steps")
        _check_required(self.qed, "qed")
        _check_required(self.rings, "rings")

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        for key in ("route_quality", "docking_std"):
            if d.get(key) is None:
                del d[key]
        return d


@dataclass
class MetricsEvent:
    molecule_count: int
    non_trivial_count: int
    trivial_count: int
    avg_binding_energy: float | None
    min_binding_energy: float | None
    avg_syn_steps: float | None
    docking_success_rate: float

    def __post_init__(self) -> None:
        if self.trivial_count + self.non_trivial_count != self.molecule_count:
            raise ValueError(
                f"trivial_count ({self.trivial_count}) + non_trivial_count "
                f"({self.non_trivial_count}) != molecule_count ({self.molecule_count})"
            )

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        for key in ("avg_binding_energy", "min_binding_energy", "avg_syn_steps"):
            if d.get(key) is None:
                del d[key]
        return d


@dataclass
class DockingProgressEvent:
    current: int
    total: int
    success_rate: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class HypothesisValidationEvent:
    hypothesis_id: str
    success: bool
    conclusion: str
    changes: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        _check_required(self.hypothesis_id, "hypothesis_id")
        _check_required(self.conclusion, "conclusion")

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        if d.get("changes") is None:
            del d["changes"]
        return d


EVENT_TYPE_MAP: dict[str, type] = {
    "start": StartEvent,
    "end": EndEvent,
    "stage": StageEvent,
    "molecule": MoleculeEvent,
    "metrics": MetricsEvent,
    "docking_progress": DockingProgressEvent,
    "hypothesis_validation": HypothesisValidationEvent,
}
```

- [ ] **Step 2: Create `__init__.py` for `src` if missing**

```bash
touch /home/dministrator/Lab/clones/molcraft-agent/src/__init__.py
```

---

### Task 2: EventLogger class

**Files:**
- Create: `src/event_logger.py`

- [ ] **Step 1: Write `src/event_logger.py`**

```python
"""EventLogger: structured JSONL logger with atomic commit and terminal passthrough.

Replaces the old StructuredLogger class. Each event type has a dedicated method
with typed parameters. The logger writes JSON Lines to a temp file and atomically
renames on commit(). Terminal output remains human-readable.
"""

import json
import os
import shutil
import sys
import tempfile
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, TextIO

from src.event_schema import (
    EndEvent,
    HypothesisValidationEvent,
    MetricsEvent,
    MoleculeEvent,
    StartEvent,
)


class EventLogger:
    """Writes structured JSONL events, passes human-readable text to terminal."""

    def __init__(self, log_path: Path) -> None:
        self.terminal: TextIO = sys.stdout
        self.log_path = log_path
        self.tmp_path = log_path.with_suffix(f".{uuid.uuid4().hex[:8]}.tmp.log")
        self.log_file = open(self.tmp_path, "w", encoding="utf-8")
        self._committed = False

    def _write_event(self, event_type: str, data: dict[str, Any]) -> None:
        event = {
            "type": event_type,
            "timestamp": datetime.now().isoformat(),
            **data,
        }
        self.log_file.write(
            json.dumps(event, ensure_ascii=False, allow_nan=False) + "\n"
        )
        self.log_file.flush()

    # ── Lifecycle ──

    def commit(self) -> None:
        self.log_file.flush()
        self.log_file.close()
        shutil.move(str(self.tmp_path), str(self.log_path))
        self._committed = True

    def close(self) -> None:
        if not self.log_file.closed:
            self.log_file.close()
        if not self._committed and self.tmp_path.exists():
            self.tmp_path.unlink(missing_ok=True)

    def flush(self) -> None:
        self.terminal.flush()
        self.log_file.flush()

    # ── Terminal passthrough (human-readable) ──

    def write(self, message: str) -> None:
        self.terminal.write(message)
        if message.strip():
            self._write_event("stdout", {"content": message})

    # ── Structured events ──

    def log_start(
        self, round: int, max_minutes: int, max_iterations: int, program_file: str
    ) -> None:
        ev = StartEvent(
            round=round,
            max_minutes=max_minutes,
            max_iterations=max_iterations,
            program_file=program_file,
        )
        self._write_event("start", ev.to_dict())

    def log_end(self, status: str) -> None:
        ev = EndEvent(status=status)
        self._write_event("end", ev.to_dict())

    def log_stage(
        self,
        name: str,
        status: str,
        duration_seconds: float | None = None,
    ) -> None:
        from src.event_schema import StageEvent

        ev = StageEvent(name=name, status=status, duration_seconds=duration_seconds)
        self._write_event("stage", ev.to_dict())

    def log_docking_progress(
        self, current: int, total: int, success_rate: float
    ) -> None:
        self._write_event(
            "docking_progress",
            {"current": current, "total": total, "success_rate": success_rate},
        )

    def log_molecule(self, event: MoleculeEvent) -> None:
        self._write_event("molecule", event.to_dict())

    def log_metrics(self, event: MetricsEvent) -> None:
        self._write_event("metrics", event.to_dict())

    def log_hypothesis_validation(
        self,
        hypothesis_id: str,
        success: bool,
        conclusion: str,
        changes: dict[str, Any] | None = None,
    ) -> None:
        ev = HypothesisValidationEvent(
            hypothesis_id=hypothesis_id,
            success=success,
            conclusion=conclusion,
            changes=changes,
        )
        self._write_event("hypothesis_validation", ev.to_dict())
```

---

### Task 3: Schema validation tests (TDD)

**Files:**
- Create: `tests/test_event_schema.py`

- [ ] **Step 1: Write failing tests**

```python
"""Tests for event dataclass schemas — field requiredness, defaults, consistency."""

import pytest
from src.event_schema import (
    MoleculeEvent,
    MetricsEvent,
    StageEvent,
    HypothesisValidationEvent,
    EVENT_TYPE_MAP,
)


class TestStageEvent:
    def test_valid_status(self):
        ev = StageEvent(name="诊断", status="begun")
        assert ev.status == "begun"

    def test_invalid_status_raises(self):
        with pytest.raises(ValueError, match="status"):
            StageEvent(name="诊断", status="invalid")

    def test_to_dict_omits_none_duration(self):
        ev = StageEvent(name="诊断", status="completed", duration_seconds=12.5)
        d = ev.to_dict()
        assert d["duration_seconds"] == 12.5

        ev2 = StageEvent(name="诊断", status="completed")
        d2 = ev2.to_dict()
        assert "duration_seconds" not in d2


class TestMoleculeEvent:
    REQUIRED = {
        "mol_smiles": "Cc1ccccc1",
        "binding_energy": -8.5,
        "composite_score": 0.85,
        "syn_steps": 2,
        "sa_score": 3.2,
        "qed": 0.7,
        "rings": 2,
    }

    def test_required_fields(self):
        ev = MoleculeEvent(**self.REQUIRED)
        assert ev.mol_smiles == "Cc1ccccc1"
        assert ev.binding_energy == -8.5

    def test_missing_required_field_raises(self):
        for field in ("mol_smiles", "binding_energy", "composite_score", "syn_steps", "qed", "rings"):
            kwargs = dict(self.REQUIRED)
            kwargs[field] = None
            with pytest.raises(ValueError, match=field):
                MoleculeEvent(**kwargs)

    def test_optional_defaults(self):
        ev = MoleculeEvent(**self.REQUIRED)
        assert ev.trivial is False
        assert ev.route_quality is None
        assert ev.docking_std is None

    def test_to_dict_omits_none_optional(self):
        ev = MoleculeEvent(**self.REQUIRED)
        d = ev.to_dict()
        assert "route_quality" not in d
        assert "docking_std" not in d
        assert "trivial" in d  # False is kept (not None)

    def test_to_dict_includes_set_optionals(self):
        ev = MoleculeEvent(**self.REQUIRED, route_quality=0.9, docking_std=0.05)
        d = ev.to_dict()
        assert d["route_quality"] == 0.9
        assert d["docking_std"] == 0.05

    def test_sa_score_can_be_none(self):
        kwargs = dict(self.REQUIRED)
        kwargs["sa_score"] = None
        ev = MoleculeEvent(**kwargs)
        assert ev.sa_score is None


class TestMetricsEvent:
    def test_counts_must_sum_to_molecule_count(self):
        MetricsEvent(
            molecule_count=10,
            non_trivial_count=8,
            trivial_count=2,
            avg_binding_energy=-7.5,
            min_binding_energy=-9.0,
            avg_syn_steps=1.5,
            docking_success_rate=0.9,
        )

    def test_count_mismatch_raises(self):
        with pytest.raises(ValueError, match="trivial_count"):
            MetricsEvent(
                molecule_count=10,
                non_trivial_count=9,
                trivial_count=2,
                avg_binding_energy=-7.5,
                min_binding_energy=-9.0,
                avg_syn_steps=1.5,
                docking_success_rate=0.9,
            )

    def test_to_dict_omits_none(self):
        ev = MetricsEvent(
            molecule_count=10,
            non_trivial_count=8,
            trivial_count=2,
            avg_binding_energy=None,
            min_binding_energy=-9.0,
            avg_syn_steps=None,
            docking_success_rate=1.0,
        )
        d = ev.to_dict()
        assert "avg_binding_energy" not in d
        assert "avg_syn_steps" not in d
        assert d["min_binding_energy"] == -9.0


class TestHypothesisValidationEvent:
    def test_required_fields(self):
        ev = HypothesisValidationEvent(hypothesis_id="H001", success=True, conclusion="ACCEPTED")
        assert ev.hypothesis_id == "H001"

    def test_missing_hypothesis_id_raises(self):
        with pytest.raises(ValueError, match="hypothesis_id"):
            HypothesisValidationEvent(hypothesis_id=None, success=True, conclusion="ACCEPTED")

    def test_changes_optional(self):
        ev = HypothesisValidationEvent(hypothesis_id="H001", success=True, conclusion="OK")
        assert ev.changes is None

        ev2 = HypothesisValidationEvent(
            hypothesis_id="H001",
            success=True,
            conclusion="OK",
            changes={"threshold": {"from": 8, "to": 6}},
        )
        assert ev2.changes == {"threshold": {"from": 8, "to": 6}}


class TestEventTypeMap:
    def test_all_events_have_entry(self):
        assert "start" in EVENT_TYPE_MAP
        assert "end" in EVENT_TYPE_MAP
        assert "stage" in EVENT_TYPE_MAP
        assert "molecule" in EVENT_TYPE_MAP
        assert "metrics" in EVENT_TYPE_MAP
        assert "docking_progress" in EVENT_TYPE_MAP
        assert "hypothesis_validation" in EVENT_TYPE_MAP
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -m pytest tests/test_event_schema.py -v 2>&1 | head -30
```
Expected: ImportError (module not found, files not created yet)

- [ ] **Step 3: Create `src/event_schema.py`** (from Task 1) and ensure the module is importable

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -c "from src.event_schema import MoleculeEvent; print('OK')"
```
Expected: OK

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -m pytest tests/test_event_schema.py -v
```
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && git add src/event_schema.py tests/test_event_schema.py src/__init__.py && git commit -m "feat: event dataclass schemas with validation"
```

---

### Task 4: EventLogger unit tests (TDD)

**Files:**
- Create: `tests/test_event_logger.py`

- [ ] **Step 1: Write failing tests**

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -m pytest tests/test_event_logger.py -v 2>&1 | head -10
```
Expected: ImportError (EventLogger not yet created)

- [ ] **Step 3: Create `src/event_logger.py`** (from Task 2)

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -c "from src.event_logger import EventLogger; print('OK')"
```
Expected: OK

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -m pytest tests/test_event_logger.py -v
```
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && git add src/event_logger.py tests/test_event_logger.py && git commit -m "feat: EventLogger with unified event type methods"
```

---

### Task 5: Integrate EventLogger into main.py

**Files:**
- Modify: `main.py`

This task replaces `StructuredLogger` class with `EventLogger`, removes the pipeline log merging block, and updates the metrics event to use the unified shape.

- [ ] **Step 1: Read current main.py to understand the full file**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && wc -l main.py
```

- [ ] **Step 2: In main.py, replace class `StructuredLogger` with import of `EventLogger`**

Replace lines that define the `StructuredLogger` class (lines 118-173) with:
```python
import sys
sys.stdout = old_stdout  # ensure we have the reference

# EventLogger is imported and used instead of StructuredLogger
from src.event_logger import EventLogger
from src.event_schema import MetricsEvent
```

- [ ] **Step 3: Replace `StructuredLogger` instantiation** with `EventLogger`

Find:
```python
tee = StructuredLogger(LOG_PATH)
```
Replace with:
```python
tee = EventLogger(LOG_PATH)
```

- [ ] **Step 4: Replace the metrics computation in main.py's finally block** (lines ~390-406)

Old code:
```python
trivial_count = sum(
    1 for m in molecules
    if m.get('route', '') == m.get('mol_smiles', '') + '>>' + m.get('mol_smiles', '')
)
trivial_ratio = trivial_count / max(len(molecules), 1)
tee.log_event("metrics",
             sample_count=len(molecules),
             trivial_ratio=trivial_ratio)
```

Replace with:
```python
trivial_count = sum(
    1 for m in molecules
    if m.get('route', '') == m.get('mol_smiles', '') + '>>' + m.get('mol_smiles', '')
)
non_trivial_count = len(molecules) - trivial_count
tee.log_metrics(MetricsEvent(
    molecule_count=len(molecules),
    non_trivial_count=non_trivial_count,
    trivial_count=trivial_count,
    avg_binding_energy=None,   # pipeline already wrote the detailed metrics
    min_binding_energy=None,
    avg_syn_steps=None,
    docking_success_rate=0.0,
))
```

- [ ] **Step 5: Remove the pipeline log merging block** (lines ~413-424)

Delete the block that reads and merges `MOLCRAFT_LOG_PATH`, since pipeline now receives EventLogger directly.

- [ ] **Step 6: Check `sys.stdout` is set to EventLogger** as in original StructuredLogger pattern

Ensure the line that redirects stdout to EventLogger instance is present:
```python
sys.stdout = tee
```

- [ ] **Step 7: Run a dry syntax check**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -c "import main; print('syntax OK')" 2>&1
```
Expected: syntax OK (may fail on other imports — ensure at least no SyntaxError)

- [ ] **Step 8: Run existing EventLogger tests to confirm no regression**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -m pytest tests/test_event_logger.py tests/test_event_schema.py -v
```
Expected: All PASS

- [ ] **Step 9: Commit**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && git add main.py && git commit -m "refactor: replace StructuredLogger with EventLogger in main.py"
```

---

### Task 6: Update pipeline.py to use EventLogger

**Files:**
- Modify: `tools/pipeline.py`

- [ ] **Step 1: Read `run_evolutionary_pipeline` signature**

Check current function signature and where `write_result_log` is called.

- [ ] **Step 2: Add `logger: EventLogger | None` parameter to `run_evolutionary_pipeline`**

```python
def run_evolutionary_pipeline(
    n_generate=50,
    n_top=10,
    strategy="mutate",
    n_generations=2,
    n_offspring_per_seed=3,
    output_dir="output",
    use_docking_guidance=True,
    append_result_log=True,
    logger=None,  # EventLogger instance (optional, for structured logging)
):
```

- [ ] **Step 3: In the docking loop, replace `sys.stdout.log_event` with `logger.log_docking_progress`**

Find calls like `sys.stdout.log_event("docking_progress"` — but these happen in `src/evaluator.py`, not pipeline.py. In pipeline.py, the docking progress events need to be emitted from `batch_dock` calls. Check the current codebase for where `docking_progress` is emitted.

The existing code likely emits docking_progress via a callback. If not, add it after the batch_dock call.

- [ ] **Step 4: In the scored_candidates loop, add molecule events**

After line ~256 (the `scored_candidates.append` block), add:
```python
if logger is not None:
    logger.log_molecule(MoleculeEvent(
        mol_smiles=smiles,
        binding_energy=mol.get("binding_energy") or 0.0,
        composite_score=candidate_composite,
        syn_steps=syn.get("steps", 0),
        sa_score=syn.get("sa_score"),
        qed=mol.get("qed") or 0.0,
        rings=count_rings(smiles),
        trivial=is_trivial,
        route_quality=quality,
        docking_std=mol.get("consensus_std"),
    ))
```

Note: `count_rings` needs to be implemented or imported from RDKit:
```python
from rdkit import Chem

def count_rings(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return Chem.rdMolDescriptors.CalcNumRings(mol)
```

- [ ] **Step 5: Replace `write_result_log` call at end of pipeline**

Find where `write_result_log` is called and replace with:
```python
if logger is not None:
    trivial_count = sum(1 for r in results if r.get("trivial", False))
    energies = [r["binding_energy"] for r in results if r.get("binding_energy") is not None]
    steps = [r.get("syn_steps", 0) for r in results]
    logger.log_metrics(MetricsEvent(
        molecule_count=len(results),
        non_trivial_count=len(results) - trivial_count,
        trivial_count=trivial_count,
        avg_binding_energy=sum(energies) / len(energies) if energies else None,
        min_binding_energy=min(energies) if energies else None,
        avg_syn_steps=sum(steps) / len(steps) if steps else None,
        docking_success_rate=0.0,
    ))

# Still write CSV (unchanged logic)
```

- [ ] **Step 6: Remove `write_result_log` function definition** (lines 43-92)

Delete the `def write_result_log(...)` function entirely.

- [ ] **Step 7: Run syntax check**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -c "from tools.pipeline import run_evolutionary_pipeline; print('syntax OK')"
```
Expected: syntax OK

- [ ] **Step 8: Run all tests**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -m pytest tests/ -v
```
Expected: All tests PASS

- [ ] **Step 9: Commit**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && git add tools/pipeline.py && git commit -m "refactor: pipeline uses EventLogger, remove write_result_log"
```

---

### Task 7: Add begin_stage/end_stage tools

**Files:**
- Create: references to tools in `tools/tools.py` (verify the correct file)

- [ ] **Step 1: Find the tools file**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && find tools -name "*.py" -maxdepth 1 | sort
```

- [ ] **Step 2: Add begin_stage/end_stage functions**

```python
import time

_stage_name: str | None = None
_stage_start: float | None = None


def begin_stage(name: str) -> str:
    """Agent tool: mark the beginning of a phase (诊断/代码演进/实验验证/复盘).

    Call this at the start of each major phase. Must be paired with end_stage().
    Nested stages are not allowed.
    """
    global _stage_name, _stage_start
    if _stage_name is not None:
        return f"错误: 阶段「{_stage_name}」尚未结束，不能开始新阶段「{name}」"
    _stage_name = name
    _stage_start = time.monotonic()
    try:
        sys.stdout.log_stage(name=name, status="begun")
    except AttributeError:
        pass  # logger not available (direct terminal run)
    return f"阶段「{name}」开始"


def end_stage() -> str:
    """Agent tool: mark the end of the current phase.

    Automatically calculates and records the duration.
    """
    global _stage_name, _stage_start
    if _stage_name is None:
        return "错误: 没有正在进行的阶段"
    elapsed = time.monotonic() - _stage_start
    name = _stage_name
    _stage_name = None
    _stage_start = None
    try:
        sys.stdout.log_stage(name=name, status="completed", duration_seconds=round(elapsed, 1))
    except AttributeError:
        pass
    return f"阶段「{name}」完成，耗时 {elapsed:.1f}s"
```

- [ ] **Step 3: Register these as tools available to the Agent**

Check how existing tools are registered (likely as functions in tools/ that the Agent calls). Ensure `begin_stage` and `end_stage` are discoverable.

- [ ] **Step 4: Run import check**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -c "from tools.tools import begin_stage, end_stage; print('OK')"
```
Expected: OK

- [ ] **Step 5: Commit**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && git add tools/tools.py && git commit -m "feat: add begin_stage/end_stage Agent tools"
```

---

### Task 8: Update program.md with stage instructions

**Files:**
- Modify: `program.md`

- [ ] **Step 1: Read program.md to find phase transition points**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && grep -n "^##" program.md
```

- [ ] **Step 2: Add `begin_stage` call at the start of each phase section**

Example additions:
```markdown
## 1. 文献分析

!!! begin_stage("文献分析")
...
```

(Exact syntax depends on how program.md instructs tool calls. Read the current format and follow the same pattern.)

- [ ] **Step 3: Add `end_stage` call at the end of each phase section**

- [ ] **Step 4: Verify no syntax errors in program.md**

program.md is markdown (not executable), so no automated check — just verify the tool call patterns match what tools.py expects.

- [ ] **Step 5: Commit**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && git add program.md && git commit -m "feat: instruct Agent to call begin_stage/end_stage at phase transitions"
```

---

### Task 9: Integration validation

- [ ] **Step 1: Run full test suite**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -m pytest tests/ -v
```
Expected: All tests PASS

- [ ] **Step 2: Dry-run the pipeline with EventLogger**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && python -c "
from pathlib import Path
from src.event_logger import EventLogger
from tools.pipeline import run_evolutionary_pipeline
import tempfile, os

tmpdir = tempfile.mkdtemp()
log_path = Path(tmpdir) / 'result.log'
logger = EventLogger(log_path)
logger.log_start(round=1, max_minutes=5, max_iterations=1, program_file='test')
logger.log_stage(name='实验验证', status='begun')

# Run pipeline with small params
try:
    run_evolutionary_pipeline(
        n_generate=5, n_top=2, n_generations=1, output_dir=tmpdir,
        logger=logger,
    )
except SystemExit:
    pass

logger.log_stage(name='实验验证', status='completed', duration_seconds=1.0)
logger.log_end(status='success')
logger.commit()

# Verify result.log has all expected event types
import json
with open(log_path) as f:
    events = [json.loads(l) for l in f if l.strip()]
types = [e['type'] for e in events]
print(f'Event types found: {set(types)}')
assert 'start' in types, 'Missing start event'
assert 'end' in types, 'Missing end event'
assert 'molecule' in types, 'Missing molecule event'
assert 'metrics' in types, 'Missing metrics event'
print('All required event types present. Integration OK.')
"
```
Expected: Prints "All required event types present. Integration OK."

- [ ] **Step 3: Final commit**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && git add -A && git status
```
Verify only expected files are modified, then commit if there are uncommitted changes.

```bash
cd /home/dministrator/Lab/clones/molcraft-agent && git commit -m "feat: integrate EventLogger across main.py, pipeline, and tools"
```
