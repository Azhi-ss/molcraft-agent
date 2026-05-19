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
