import pytest
from src.calibrator import fit_binding_score, fit_sa_score, CalibrationResult


class TestFitBindingScore:
    def test_perfect_clipped_linear(self):
        # vina=-10 → 1.0, vina=-5 → 0.5, vina=0 → 0.0
        # implies threshold≈0, range≈10
        points = [
            {"vina_raw": -10.0, "actual_score": 1.0},
            {"vina_raw": -5.0, "actual_score": 0.5},
            {"vina_raw": 0.0, "actual_score": 0.0},
        ]
        result = fit_binding_score(points)
        assert result.function == "clipped_linear"
        assert result.params["threshold"] == pytest.approx(0.0, abs=0.5)
        assert result.params["range"] == pytest.approx(10.0, abs=1.0)

    def test_single_point(self):
        points = [{"vina_raw": -9.591, "actual_score": 0.1766}]
        result = fit_binding_score(points)
        assert result.function in ("clipped_linear", "minmax")
        # Verify the point maps correctly
        from src.scorer import compute_binding_score
        cal = {"binding_score": {"function": result.function, "params": result.params}}
        predicted = compute_binding_score(-9.591, cal)
        assert predicted == pytest.approx(0.1766, abs=0.05)

    def test_empty_returns_defaults(self):
        result = fit_binding_score([])
        assert result.params["threshold"] == 0.0
        assert result.params["range"] == 15.0


class TestFitSAScore:
    def test_from_known_data(self):
        # sa=2.0 → 0.5, sa=3.0 → 0.25, sa=4.0 → 0.0
        # implies cutoff≈4, scale≈4
        points = [
            {"sa_raw": 2.0, "actual_score": 0.5},
            {"sa_raw": 3.0, "actual_score": 0.25},
            {"sa_raw": 4.0, "actual_score": 0.0},
        ]
        result = fit_sa_score(points)
        assert result.params["cutoff"] == pytest.approx(4.0, abs=0.5)
        assert result.params["scale"] == pytest.approx(4.0, abs=0.5)

    def test_empty_returns_defaults(self):
        result = fit_sa_score([])
        assert result.params["cutoff"] == 4.0
        assert result.params["scale"] == 4.0
