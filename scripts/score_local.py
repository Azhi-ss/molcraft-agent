#!/usr/bin/env python3
"""Local scoring CLI — score a result.csv or calibrate normalization params.

Usage:
    python scripts/score_local.py output/result.csv --vina output/result.json
    python scripts/score_local.py --calibrate --vina output/result.json --actual binding_score=0.1766 sa_score=0.656
"""
import argparse
import json
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.scorer import score_csv, load_calibration, save_calibration
from src.calibrator import calibrate_from_submission


def load_vina_scores_from_json(path: str) -> dict[str, float]:
    """Load Vina scores from result.json or molecules.jsonl.

    Returns dict mapping SMILES -> binding energy (kcal/mol).
    """
    with open(path) as f:
        data = json.load(f)

    scores = {}
    # result.json with "molecules" key
    if isinstance(data, dict) and "molecules" in data:
        mols = data["molecules"]
        if isinstance(mols, list):
            for mol in mols:
                if isinstance(mol, dict) and "smiles" in mol and "be" in mol:
                    scores[mol["smiles"]] = mol["be"]
    # result.json with "top_molecules" key
    if isinstance(data, dict) and "top_molecules" in data:
        for mol in data["top_molecules"]:
            if isinstance(mol, dict) and "smiles" in mol and "be" in mol:
                scores[mol["smiles"]] = mol["be"]

    # Fallback: try molecules.jsonl in same directory
    if not scores:
        jsonl_path = os.path.join(os.path.dirname(path), "molecules.jsonl")
        if os.path.exists(jsonl_path):
            with open(jsonl_path) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    entry = json.loads(line)
                    if "molecules" in entry:
                        for mol in entry["molecules"]:
                            if "smiles" in mol and "be" in mol:
                                scores[mol["smiles"]] = mol["be"]

    return scores


def main():
    parser = argparse.ArgumentParser(description="Local scoring system for MolCraft competition")
    parser.add_argument("csv_path", nargs="?", help="Path to result.csv")
    parser.add_argument("--vina", help="Path to result.json or molecules.jsonl for Vina scores")
    parser.add_argument("--calibrate", action="store_true", help="Calibrate normalization params from submission data")
    parser.add_argument("--actual", nargs="*", help="Actual scores from submission: key=value pairs (e.g. binding_score=0.1766)")
    parser.add_argument("--calibration", help="Path to calibration.json (default: data/calibration.json)")

    args = parser.parse_args()

    if args.calibrate:
        if not args.vina or not args.actual:
            print("Error: --calibrate requires --vina and --actual")
            sys.exit(1)

        vina_scores = load_vina_scores_from_json(args.vina)
        if not vina_scores:
            print(f"Error: No Vina scores found in {args.vina}")
            sys.exit(1)

        actual_scores = {}
        for pair in args.actual:
            k, v = pair.split("=")
            actual_scores[k] = float(v)

        from rdkit import Chem
        from src.evaluator import estimate_sa_score
        sa_raws = {}
        for smi in vina_scores:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                sa_raws[smi] = estimate_sa_score(mol)

        cal = calibrate_from_submission(vina_scores, sa_raws, actual_scores)
        save_calibration(cal, args.calibration)
        print("Calibration updated and saved to data/calibration.json:")
        print(json.dumps(cal, indent=2))
        return

    if not args.csv_path:
        parser.print_help()
        sys.exit(1)

    # Scoring mode
    vina_scores = {}
    if args.vina:
        vina_scores = load_vina_scores_from_json(args.vina)

    cal = load_calibration(args.calibration)
    result = score_csv(args.csv_path, vina_scores, cal)

    print("=" * 50)
    print("Local Score Report")
    print("=" * 50)
    for key, value in sorted(result.items()):
        if isinstance(value, float):
            print(f"  {key:40s}: {value:.6f}")
        else:
            print(f"  {key:40s}: {value}")
    print("=" * 50)


if __name__ == "__main__":
    main()