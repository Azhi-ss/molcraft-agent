"""Event type dataclasses for structured JSONL logging.

Each event type is a dataclass with a to_dict() method for JSON serialization.
Validation occurs at construction time.
"""

from dataclasses import dataclass, asdict
from typing import Any


def _check_required(value: Any, field_name: str) -> None:
    if value is None:
        raise ValueError(f"{field_name} is required")
    if isinstance(value, str) and not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string")


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
