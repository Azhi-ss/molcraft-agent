"""Tests for scorer.py - validity, SA normalization, and binding normalization."""

import pytest
from pathlib import Path

# Import the functions to test
from src.scorer import (
    load_calibration,
    save_calibration,
    compute_validity_score,
    compute_sa_score_normalized,
    compute_binding_score,
    compute_route_validity_score,
    compute_balance_score,
    compute_step_penalty_score,
    compute_starting_material_availability_score,
    compute_mol_score,
    compute_route_score,
    compute_total_score,
    _is_trivial_route,
    score_csv,
)


class TestValidityScore:
    """Test SMILES validity scoring."""

    def test_valid_molecule(self):
        """Valid SMILES returns 1.0."""
        assert compute_validity_score("CCO") == 1.0

    def test_invalid_smiles(self):
        """Invalid SMILES returns 0.0."""
        assert compute_validity_score("InvalidSMILES") == 0.0

    def test_empty_string(self):
        """Empty string returns 0.0."""
        assert compute_validity_score("") == 0.0


class TestSAScoreNormalized:
    """Test SA score normalization."""

    def test_sa_above_cutoff(self):
        """SA >= cutoff returns 0.0."""
        assert compute_sa_score_normalized(5.0) == 0.0
        assert compute_sa_score_normalized(4.0) == 0.0

    def test_sa_at_cutoff(self):
        """SA at cutoff returns 0.0."""
        assert compute_sa_score_normalized(4.0) == 0.0

    def test_sa_midpoint(self):
        """SA at midpoint (2.0) returns 0.5."""
        assert compute_sa_score_normalized(2.0) == 0.5

    def test_sa_zero(self):
        """SA = 0 returns 1.0."""
        assert compute_sa_score_normalized(0.0) == 1.0

    def test_sa_negative_clamped(self):
        """SA < 0 is clamped to 1.0."""
        assert compute_sa_score_normalized(-1.0) == 1.0


class TestBindingScore:
    """Test binding score (Vina) normalization."""

    def test_vina_negative_nine(self):
        """vina=-9, threshold=0, range=15: (0-(-9))/15 = 0.6."""
        assert compute_binding_score(-9.0) == 0.6

    def test_vina_zero(self):
        """vina=0 returns 0.0."""
        assert compute_binding_score(0.0) == 0.0

    def test_vina_minus_fifteen(self):
        """vina=-15: (0-(-15))/15 = 1.0."""
        assert compute_binding_score(-15.0) == 1.0

    def test_vina_positive_clamped(self):
        """vina > threshold is clamped to 0.0."""
        assert compute_binding_score(5.0) == 0.0


class TestRouteValidityScore:
    def test_all_valid(self):
        assert compute_route_validity_score([True] * 10) == 1.0

    def test_half_valid(self):
        assert compute_route_validity_score([True] * 5 + [False] * 5) == 0.5

    def test_empty(self):
        assert compute_route_validity_score([]) == 0.0


class TestBalanceScore:
    def test_perfect_balance(self):
        assert compute_balance_score(12, 12) == 1.0

    def test_one_off(self):
        score = compute_balance_score(12, 11)
        assert 0.8 < score <= 1.0

    def test_large_imbalance(self):
        # ratio=8/12=0.667 > 0.5 so not zero; 12->4 (ratio=0.333 < 0.5) triggers 0.0
        assert compute_balance_score(12, 4) == 0.0


class TestStepPenaltyScore:
    def test_one_step(self):
        assert compute_step_penalty_score(1) == 1.0

    def test_three_steps(self):
        assert compute_step_penalty_score(3) == pytest.approx(0.7)

    def test_five_steps(self):
        # 1.0 - 0.15*(5-1) = 0.4
        assert compute_step_penalty_score(5) == pytest.approx(0.4)


class TestMolScore:
    def test_composition(self):
        # binding=0.5, validity=1.0, sa=0.25 → 0.8*0.5 + 0.1*1.0 + 0.1*0.25 = 0.525
        assert compute_mol_score(0.5, 1.0, 0.25) == pytest.approx(0.525)


class TestRouteScore:
    def test_composition(self):
        # route_validity=1.0, sm_avail=0.9, step_penalty=0.8, convergence=1.0, balance=1.0
        # = 0.55 + 0.27 + 0.04 + 0.05 + 0.05 = 0.96
        assert compute_route_score(1.0, 0.9, 0.8, 1.0, 1.0) == pytest.approx(0.96)


class TestTotalScore:
    def test_composition(self):
        # mol=0.525, route=0.96 → 0.7*0.525 + 0.3*0.96 = 0.6555
        assert compute_total_score(0.525, 0.96) == pytest.approx(0.6555)


class TestIsTrivialRoute:
    def test_trivial_same_molecule(self):
        assert _is_trivial_route("c1ccccc1>>c1ccccc1", "c1ccccc1") == True

    def test_nontrivial(self):
        assert (
            _is_trivial_route(
                "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1", "c1ccccc1-c1ccccc1"
            )
            == False
        )

    def test_no_arrow(self):
        assert _is_trivial_route("c1ccccc1", "c1ccccc1") == True


class TestScoreCSV:
    def test_basic_scoring(self, tmp_path):
        # Write a minimal CSV, score it
        csv_file = tmp_path / "result.csv"
        csv_file.write_text(
            "mol_smiles,route\nc1ccccc1,c1ccccc1>>c1ccccc1\nc1ccc(-c2ccccc2)cc1,Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1\n"
        )
        vina = {"c1ccccc1": -5.0, "c1ccc(-c2ccccc2)cc1": -9.0}
        result = score_csv(str(csv_file), vina_scores=vina)
        assert "total_score" in result
        assert "mol_score" in result
        assert "route_score" in result
        assert "binding_score" in result
        assert "route_validity_score" in result
        assert 0.0 <= result["total_score"] <= 1.0
        assert result["sample_count"] == 2

    def test_missing_vina_gets_zero_binding(self, tmp_path):
        csv_file = tmp_path / "result.csv"
        csv_file.write_text("mol_smiles,route\nc1ccccc1,c1ccccc1>>c1ccccc1\n")
        result = score_csv(str(csv_file), vina_scores={})
        assert result["binding_score"] == 0.0
