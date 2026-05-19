"""Scorer core functions: validity, SA normalization, binding normalization."""

import json
from pathlib import Path
from rdkit import Chem


# Path to calibration file, relative to project root
CALIBRATION_PATH = Path(__file__).parent.parent / "data" / "calibration.json"


# Default calibration values
DEFAULT_CALIBRATION = {
    "binding_score": {
        "type": "clipped_linear",
        "threshold": 0.0,
        "range": 15.0,
    },
    "sa_score": {
        "type": "step",
        "cutoff": 4.0,
        "scale": 4.0,
    },
}


def load_calibration(path=None) -> dict:
    """Load calibration from JSON file.

    Args:
        path: Path to calibration JSON. Defaults to data/calibration.json.

    Returns:
        Calibration dict with binding_score and sa_score keys.
    """
    if path is None:
        path = CALIBRATION_PATH

    if Path(path).exists():
        with open(path) as f:
            return json.load(f)
    return DEFAULT_CALIBRATION.copy()


def save_calibration(cal: dict, path=None):
    """Save calibration to JSON file.

    Args:
        cal: Calibration dict to save.
        path: Path to write. Defaults to data/calibration.json.
    """
    if path is None:
        path = CALIBRATION_PATH

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(cal, f, indent=2)


def compute_validity_score(smiles: str) -> float:
    """Compute validity score for a SMILES string.

    Args:
        smiles: SMILES string to validate.

    Returns:
        1.0 if valid SMILES, 0.0 otherwise.
    """
    if not smiles:
        return 0.0
    mol = Chem.MolFromSmiles(smiles)
    return 1.0 if mol is not None else 0.0


def compute_sa_score_normalized(sa_raw: float, cal: dict = None) -> float:
    """Normalize SA score using step function.

    Args:
        sa_raw: Raw SA score (0-10, lower=better).
        cal: Optional calibration dict. Defaults to loaded calibration.

    Returns:
        Normalized score: 0.0 if sa_raw >= cutoff, else (cutoff - sa_raw) / scale,
        clamped to [0, 1].
    """
    if cal is None:
        cal = load_calibration()

    sa_cal = cal.get("sa_score", DEFAULT_CALIBRATION["sa_score"])
    cutoff = sa_cal.get("cutoff", 4.0)
    scale = sa_cal.get("scale", 4.0)

    if sa_raw >= cutoff:
        return 0.0

    score = (cutoff - sa_raw) / scale
    return max(0.0, min(1.0, score))


def compute_binding_score(vina_raw: float, cal: dict = None) -> float:
    """Normalize binding score (Vina) using clipped linear or minmax.

    Args:
        vina_raw: Raw Vina binding score (typically negative).
        cal: Optional calibration dict. Defaults to loaded calibration.

    Returns:
        Normalized score in [0, 1].
    """
    if cal is None:
        cal = load_calibration()

    binding_cal = cal.get("binding_score", DEFAULT_CALIBRATION["binding_score"])
    norm_type = binding_cal.get("type", "clipped_linear")

    if norm_type == "clipped_linear":
        threshold = binding_cal.get("threshold", 0.0)
        range_ = binding_cal.get("range", 15.0)
        score = (threshold - vina_raw) / range_
    elif norm_type == "minmax":
        min_val = binding_cal.get("min", -15.0)
        max_val = binding_cal.get("max", 0.0)
        score = (vina_raw - min_val) / (max_val - min_val)
    else:
        raise ValueError(f"Unknown binding normalization type: {norm_type}")

    return max(0.0, min(1.0, score))
