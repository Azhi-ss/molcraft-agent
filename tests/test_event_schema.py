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
