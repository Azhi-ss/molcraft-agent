"""Scorer core functions: validity, SA normalization, binding normalization."""

import csv
import json
from pathlib import Path

from rdkit import Chem

from src.evaluator import estimate_sa_score


# Path to calibration file, relative to project root
CALIBRATION_PATH = Path(__file__).parent.parent / "data" / "calibration.json"


# Default calibration values
DEFAULT_CALIBRATION = {
    "version": 1,
    "binding_score": {
        "function": "clipped_linear",
        "params": {"threshold": 0.0, "range": 15.0},
    },
    "sa_score": {
        "function": "step",
        "params": {"cutoff": 4.0, "scale": 4.0},
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
            cal = json.load(f)
        # Merge with defaults for any missing keys
        result = DEFAULT_CALIBRATION.copy()
        result.update(cal)
        for key in DEFAULT_CALIBRATION:
            if key not in cal and isinstance(DEFAULT_CALIBRATION[key], dict):
                result[key] = DEFAULT_CALIBRATION[key]
            elif key in cal and isinstance(DEFAULT_CALIBRATION.get(key), dict) and isinstance(cal[key], dict):
                merged = DEFAULT_CALIBRATION[key].copy()
                merged.update(cal[key])
                result[key] = merged
        return result
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
    """Normalize SA score to 0-1 (lower SA = higher score).

    Supports:
    - step: SA >= cutoff → 0, else (cutoff - sa) / scale
    - inverted: (max_sa - sa) / scale
    """
    if cal is None:
        cal = load_calibration()

    sa_cal = cal.get("sa_score", DEFAULT_CALIBRATION["sa_score"])
    func = sa_cal.get("function", "step")
    params = sa_cal.get("params", sa_cal)

    if func == "step":
        cutoff = params.get("cutoff", 4.0)
        scale = params.get("scale", 4.0)
        if sa_raw >= cutoff:
            return 0.0
        score = (cutoff - sa_raw) / scale
    elif func == "inverted":
        max_sa = params.get("max_sa", 10.0)
        scale = params.get("scale", 10.0)
        score = (max_sa - sa_raw) / scale
    else:
        raise ValueError(f"Unknown sa_score function: {func}")

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
    func = binding_cal.get("function", binding_cal.get("type", "clipped_linear"))
    params = binding_cal.get("params", binding_cal)  # support both nested and flat

    if func == "clipped_linear":
        threshold = params.get("threshold", 0.0)
        range_ = params.get("range", 15.0)
        score = (threshold - vina_raw) / range_
    elif func == "minmax":
        min_val = params.get("min", -15.0)
        max_val = params.get("max", 0.0)
        score = (vina_raw - min_val) / (max_val - min_val)
    else:
        raise ValueError(f"Unknown binding normalization function: {func}")

    return max(0.0, min(1.0, score))


def compute_route_validity_score(routes_valid: list[bool]) -> float:
    """Fraction of routes that pass validity checks. Empty list → 0.0."""
    if not routes_valid:
        return 0.0
    return sum(routes_valid) / len(routes_valid)


def compute_balance_score(reactant_heavy: int, product_heavy: int) -> float:
    """0-1 score based on atom balance. Perfect=1.0, >50% imbalance=0.0.
    ratio = product_heavy / reactant_heavy
    If ratio < 0.5 or ratio > 1.5 → 0.0
    Else: max(0.0, 1.0 - 2.0 * abs(1.0 - ratio))
    """
    if reactant_heavy == 0:
        return 0.0
    ratio = product_heavy / reactant_heavy
    if ratio < 0.5 or ratio > 1.5:
        return 0.0
    return max(0.0, 1.0 - 2.0 * abs(1.0 - ratio))


def compute_step_penalty_score(n_steps: int) -> float:
    """1.0 for 1 step, decreasing by 0.15 per extra step. n_steps<=0 → 0.0."""
    if n_steps <= 0:
        return 0.0
    return max(0.0, 1.0 - 0.15 * (n_steps - 1))


def compute_starting_material_availability_score(reactant_smiles: list[str]) -> float:
    """Heuristic: SA < 4 → 1.0, SA 4-6 → 0.5, SA > 6 → 0.0. Average over all reactants.
    Empty list → 0.0. Uses estimate_sa_score from evaluator.
    """
    if not reactant_smiles:
        return 0.0
    scores = []
    for smi in reactant_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        sa = estimate_sa_score(mol)
        if sa < 4:
            scores.append(1.0)
        elif sa <= 6:
            scores.append(0.5)
        else:
            scores.append(0.0)
    if not scores:
        return 0.0
    return sum(scores) / len(scores)


def compute_mol_score(binding_score: float, validity_score: float, sa_score: float) -> float:
    """Molecular score = 0.8*binding + 0.1*validity + 0.1*sa."""
    return 0.8 * binding_score + 0.1 * validity_score + 0.1 * sa_score


def compute_route_score(
    route_validity: float,
    sm_availability: float,
    step_penalty: float,
    convergence: float,
    balance: float,
) -> float:
    """Route score = 0.55*route_validity + 0.30*sm_availability + 0.05*step_penalty + 0.05*convergence + 0.05*balance."""
    return (
        0.55 * route_validity
        + 0.30 * sm_availability
        + 0.05 * step_penalty
        + 0.05 * convergence
        + 0.05 * balance
    )


def compute_total_score(mol_score: float, route_score: float) -> float:
    """Total = 0.7*mol + 0.3*route."""
    return 0.7 * mol_score + 0.3 * route_score


def _is_trivial_route(route: str, target_smiles: str) -> bool:
    """Check if route is trivial (A>>A).

    - No >> in route → True
    - Single reactant identical to product → True
    - Single reactant identical to first product (before any .) → True
    - Otherwise → False
    """
    if ">>" not in route:
        return True

    parts = route.split(">>")
    if len(parts) != 2:
        return False

    reactants_str, products_str = parts

    # Check single reactant identical to product
    reactants = [r.strip() for r in reactants_str.split(".")]
    if len(reactants) == 1:
        reactant = reactants[0]
        product = products_str.strip()
        # Remove any . separated additional products for comparison
        first_product = product.split(".")[0].strip()
        if reactant == product or reactant == first_product:
            return True

    return False


def score_csv(csv_path: str, vina_scores: dict = None, cal: dict = None) -> dict:
    """Score a result.csv file and return full scoring report.

    Args:
        csv_path: Path to result.csv (columns: mol_smiles,route)
        vina_scores: Dict mapping SMILES → Vina raw score (kcal/mol, negative)
        cal: Calibration params (loaded from file if None)

    Returns:
        Dict with keys: total_score, mol_score, route_score, binding_score,
        validity_score, sa_score, route_validity_score,
        starting_material_availability_score, sample_count
        All float values rounded to 6 decimal places.
    """
    import math

    if vina_scores is None:
        vina_scores = {}
    if cal is None:
        cal = load_calibration()

    rows = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    if not rows:
        return {
            "total_score": 0.0,
            "mol_score": 0.0,
            "route_score": 0.0,
            "binding_score": 0.0,
            "validity_score": 0.0,
            "sa_score": 0.0,
            "route_validity_score": 0.0,
            "starting_material_availability_score": 0.0,
            "sample_count": 0,
        }

    validity_scores = []
    binding_scores = []
    sa_scores = []
    route_validity_flags = []
    sm_availabilities = []
    step_penalties = []

    for row in rows:
        mol_smiles = row.get("mol_smiles", "").strip()
        route = row.get("route", "").strip()

        # Validity
        val_score = compute_validity_score(mol_smiles)
        validity_scores.append(val_score)

        # Binding
        vina_raw = vina_scores.get(mol_smiles, None)
        if vina_raw is not None:
            b_score = compute_binding_score(vina_raw, cal)
        else:
            b_score = 0.0
        binding_scores.append(b_score)

        # SA score
        mol = Chem.MolFromSmiles(mol_smiles) if mol_smiles else None
        if mol is not None:
            sa_raw = estimate_sa_score(mol)
            s_score = compute_sa_score_normalized(sa_raw, cal)
        else:
            s_score = 0.0
        sa_scores.append(s_score)

        # Route validity (trivial check)
        is_non_trivial = not _is_trivial_route(route, mol_smiles)
        route_validity_flags.append(is_non_trivial)

        # Extract reactants from LAST step for SM availability
        # Route format: step1 | step2 | step3  or  reactant1.reactant2>>product
        # Multi-step: "A>>B | C>>D | D>>E" → last step is "D>>E"
        # Single step: "A.B>>C"
        reactants_for_sm = []
        if ">>" in route:
            # Split by >> to get the last reaction step
            steps = route.split("|")
            last_step = steps[-1].strip()
            if ">>" in last_step:
                reactants_side = last_step.split(">>")[0].strip()
                reactants_for_sm = [r.strip() for r in reactants_side.split(".") if r.strip()]

        sm_avail = compute_starting_material_availability_score(reactants_for_sm)
        sm_availabilities.append(sm_avail)

        # Step count: number of | separators + 1
        n_steps = route.count("|") + 1
        step_penalty = compute_step_penalty_score(n_steps)
        step_penalties.append(step_penalty)

    # Average across molecules
    n = len(rows)
    avg_validity = sum(validity_scores) / n if n > 0 else 0.0
    avg_binding = sum(binding_scores) / n if n > 0 else 0.0
    avg_sa = sum(sa_scores) / n if n > 0 else 0.0
    avg_sm_avail = sum(sm_availabilities) / n if n > 0 else 0.0
    avg_step_penalty = sum(step_penalties) / n if n > 0 else 0.0

    # Route validity: fraction of non-trivial valid routes
    # Use route_validity_flags directly (True = non-trivial = valid route)
    route_validity = sum(route_validity_flags) / n if n > 0 else 0.0

    # Compute balance scores from actual route atom counts
    balance_scores = []
    for row in rows:
        mol_smiles = row.get("mol_smiles", "").strip()
        route = row.get("route", "").strip()
        if ">>" not in route or _is_trivial_route(route, mol_smiles):
            balance_scores.append(0.0)
            continue
        step_scores = []
        for step in route.split("|"):
            step = step.strip()
            if ">>" not in step:
                continue
            left, right = step.split(">>", 1)
            reactants = [r.strip() for r in left.split(".") if r.strip()]
            products = [p.strip() for p in right.split(".") if p.strip()]
            r_heavy = 0
            p_heavy = 0
            for r in reactants:
                m = Chem.MolFromSmiles(r)
                if m:
                    r_heavy += m.GetNumHeavyAtoms()
            for p in products:
                m = Chem.MolFromSmiles(p)
                if m:
                    p_heavy += m.GetNumHeavyAtoms()
            step_scores.append(compute_balance_score(r_heavy, p_heavy))
        if step_scores:
            balance_scores.append(sum(step_scores) / len(step_scores))
        else:
            balance_scores.append(0.0)
    avg_balance = sum(balance_scores) / len(balance_scores) if balance_scores else 0.0

    # Route validity: non-trivial AND all steps have non-zero balance
    for i, row in enumerate(rows):
        if balance_scores[i] == 0.0 and route_validity_flags[i]:
            route_validity_flags[i] = False
    route_validity = sum(route_validity_flags) / n if n > 0 else 0.0

    # Use placeholder for convergence (route convergence not easy to compute)
    convergence = 1.0

    # Compute sub-scores
    mol_score = compute_mol_score(avg_binding, avg_validity, avg_sa)
    route_score = compute_route_score(
        route_validity, avg_sm_avail, avg_step_penalty, convergence, avg_balance
    )
    total_score = compute_total_score(mol_score, route_score)

    def r6(x):
        """Round to 6 decimal places."""
        return math.floor(x * 1_000_000 + 0.5) / 1_000_000

    return {
        "total_score": r6(total_score),
        "mol_score": r6(mol_score),
        "route_score": r6(route_score),
        "binding_score": r6(avg_binding),
        "validity_score": r6(avg_validity),
        "sa_score": r6(avg_sa),
        "route_validity_score": r6(route_validity),
        "starting_material_availability_score": r6(avg_sm_avail),
        "sample_count": n,
    }
