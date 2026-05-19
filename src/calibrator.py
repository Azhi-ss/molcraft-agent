"""Fitting of normalization parameters from (raw_score, actual_score) pairs."""

from dataclasses import dataclass


@dataclass
class CalibrationResult:
    function: str
    params: dict
    residuals: list[float]
    rms_error: float


def fit_binding_score(points: list[dict], function: str = "clipped_linear") -> CalibrationResult:
    """Fit binding_score normalization from (vina_raw, actual_score) pairs.

    For clipped_linear: score = clamp((threshold - vina) / range, 0, 1)
    Linear regression: vina = threshold - score * range
    → vina = a + b*score where a=threshold, b=-range
    Use least squares on (score, vina) to find a and b.

    For minmax: score = (vina - min) / (max - min)
    Just use min/max of observed vina values.

    Special cases:
    - Single point: assume threshold=0, solve for range = (0 - vina) / score
    - All scores identical: return defaults (threshold=0, range=15)
    - Empty points: return defaults
    """
    # Default calibration
    default_params = {"threshold": 0.0, "range": 15.0}

    if not points:
        return CalibrationResult(function=function, params=default_params, residuals=[], rms_error=0.0)

    # Filter points where score is in [0, 1] — boundary points constrain threshold
    active = [(p["vina_raw"], p["actual_score"]) for p in points if 0 <= p["actual_score"] <= 1]

    if len(active) == 1:
        # Single point: assume threshold=0, solve for range = (0 - vina) / score
        vina, score = active[0]
        if score > 0:
            range_ = (0.0 - vina) / score
        else:
            range_ = 15.0
        threshold = 0.0
        params = {"threshold": round(threshold, 6), "range": round(range_, 6)}
        residuals = []
        rms_error = 0.0
        return CalibrationResult(function=function, params=params, residuals=residuals, rms_error=rms_error)

    if len(active) < 2:
        # Fall back to defaults
        return CalibrationResult(function=function, params=default_params, residuals=[], rms_error=0.0)

    # Least-squares fit: vina = a + b * score
    n = len(active)
    sum_x = sum(s for _, s in active)
    sum_y = sum(v for v, _ in active)
    sum_xy = sum(v * s for v, s in active)
    sum_x2 = sum(s * s for _, s in active)

    denom = n * sum_x2 - sum_x * sum_x
    if abs(denom) < 1e-12:
        return CalibrationResult(function=function, params=default_params, residuals=[], rms_error=0.0)

    b = (n * sum_xy - sum_x * sum_y) / denom
    a = (sum_y - b * sum_x) / n

    threshold = a
    range_ = -b

    if range_ <= 0:
        range_ = 15.0
        threshold = 0.0

    params = {"threshold": round(threshold, 6), "range": round(range_, 6)}

    # Compute residuals and RMS error
    residuals = []
    for vina, score in active:
        predicted = max(0.0, min(1.0, (params["threshold"] - vina) / params["range"]))
        residuals.append(predicted - score)
    rms_error = (sum(r * r for r in residuals) / len(residuals)) ** 0.5

    return CalibrationResult(function=function, params=params, residuals=residuals, rms_error=rms_error)


def fit_sa_score(points: list[dict]) -> CalibrationResult:
    """Fit sa_score normalization from (sa_raw, actual_score) pairs.

    Model: score = max(0, (cutoff - sa) / scale)
    - cutoff is between max nonzero sa and min zero sa
    - scale = average of (cutoff - sa) / score for nonzero points
    - Empty points: return defaults (cutoff=4, scale=4)
    """
    default_params = {"cutoff": 4.0, "scale": 4.0}

    if not points:
        return CalibrationResult(function="step", params=default_params, residuals=[], rms_error=0.0)

    # Separate nonzero and zero-score points
    nonzero = [(p["sa_raw"], p["actual_score"]) for p in points if p["actual_score"] > 0]
    zero = [p["sa_raw"] for p in points if p["actual_score"] == 0]

    if not nonzero:
        return CalibrationResult(function="step", params=default_params, residuals=[], rms_error=0.0)

    # cutoff is at max zero sa (the highest sa with zero score)
    max_zero = max(zero) if zero else None

    if max_zero is not None:
        # Use the highest zero-score sa as cutoff
        cutoff = max_zero
    else:
        # No zero-score points — use midpoint between min and max nonzero sa
        all_sa = [s for s, _ in nonzero]
        cutoff = (min(all_sa) + max(all_sa)) / 2

    # scale = average of (cutoff - sa) / score for nonzero points
    scales = [(cutoff - sa) / score for sa, score in nonzero if score > 0]
    if scales:
        scale = sum(scales) / len(scales)
    else:
        scale = 4.0

    if scale <= 0:
        scale = 4.0

    params = {"cutoff": round(cutoff, 6), "scale": round(scale, 6)}

    # Compute residuals and RMS error
    residuals = []
    for p in points:
        sa_raw = p["sa_raw"]
        actual = p["actual_score"]
        if sa_raw >= cutoff:
            predicted = 0.0
        else:
            predicted = max(0.0, min(1.0, (params["cutoff"] - sa_raw) / params["scale"]))
        residuals.append(predicted - actual)
    rms_error = (sum(r * r for r in residuals) / len(residuals)) ** 0.5

    return CalibrationResult(function="step", params=params, residuals=residuals, rms_error=rms_error)


def calibrate_from_submission(
    vina_scores: dict[str, float],
    sa_raws: dict[str, float],
    actual_scores: dict,
) -> dict:
    """Full calibration from a single submission result.

    Args:
        vina_scores: SMILES → Vina raw score
        sa_raws: SMILES → SAScore raw
        actual_scores: dict with competition scores (binding_score, sa_score, etc.)

    Returns:
        Updated calibration dict suitable for save_calibration()
    """
    # Build (vina_raw, actual_score) points from vina_scores
    vina_points = [
        {"vina_raw": v, "actual_score": actual_scores.get(smiles, 0.0)}
        for smiles, v in vina_scores.items()
    ]
    sa_points = [
        {"sa_raw": sa_raws.get(smiles, 10.0), "actual_score": actual_scores.get(smiles, 0.0)}
        for smiles in vina_scores
    ]

    binding_result = fit_binding_score(vina_points)
    sa_result = fit_sa_score(sa_points)

    return {
        "version": 1,
        "binding_score": {
            "function": binding_result.function,
            "params": binding_result.params,
        },
        "sa_score": {
            "function": sa_result.function,
            "params": sa_result.params,
        },
    }
