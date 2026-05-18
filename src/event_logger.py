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
