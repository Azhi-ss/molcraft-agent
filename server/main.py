#!/usr/bin/env python3
"""PocketXMol GPU Inference Server.

FastAPI service that exposes PocketXMol molecule generation via HTTP API.
Deploy on a GPU server and set DIFFUSION_API_URL on the Agent machine.

Usage:
    cd /root/PocketXMol
    pip install fastapi uvicorn
    python /path/to/server/main.py --pxm-dir /root/PocketXMol --port 8000

Endpoints:
    GET  /health  - Health check
    POST /generate - Generate molecules for a protein pocket
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
import yaml
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

app = FastAPI(title="PocketXMol Inference Server")

# Configured at startup
_PXM_DIR: str = ""
_DEVICE: str = "cuda:0"
_CHECKPOINT_EXISTS: bool = False


class GenerateRequest(BaseModel):
    pdb_content: str = Field(description="Protein PDB file content")
    n_molecules: int = Field(default=20, description="Number of molecules to generate")
    pocket_center: list[float] = Field(
        default=[0.0, 0.0, 0.0],
        description="Pocket center coordinates [x, y, z]",
    )
    pocket_radius: float = Field(default=15.0, description="Pocket radius in Angstroms")
    task: str = Field(default="dock_smallmol", description="Task type")


class MoleculeResult(BaseModel):
    smiles: str
    score: float = 0.0
    qed: Optional[float] = None
    mw: Optional[float] = None
    logp: Optional[float] = None
    sa_score: Optional[float] = None


class GenerateResponse(BaseModel):
    status: str
    molecules: list[MoleculeResult]
    generation_time_seconds: float
    message: str = ""


class HealthResponse(BaseModel):
    status: str
    model_loaded: bool
    gpu_available: bool
    device: str = ""
    pxm_dir: str = ""


@app.get("/health", response_model=HealthResponse)
def health_check():
    import torch
    return HealthResponse(
        status="ok" if _CHECKPOINT_EXISTS else "checkpoint_missing",
        model_loaded=_CHECKPOINT_EXISTS,
        gpu_available=torch.cuda.is_available(),
        device=_DEVICE,
        pxm_dir=_PXM_DIR,
    )


@app.post("/generate", response_model=GenerateResponse)
def generate_molecules(request: GenerateRequest):
    if not _CHECKPOINT_EXISTS:
        raise HTTPException(
            status_code=503,
            detail=f"Model checkpoint not found at {_PXM_DIR}/data/trained_models/pxm/checkpoints/pocketxmol.ckpt",
        )

    start_time = time.time()
    outdir = None

    try:
        # Write PDB to temp file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False, prefix="pxm_input_"
        ) as f:
            f.write(request.pdb_content)
            pdb_path = f.name

        # Create temp output dir
        outdir = tempfile.mkdtemp(prefix="pxm_output_")

        # Generate task-specific YAML config
        config_path = _write_task_config(
            pdb_path=pdb_path,
            outdir=outdir,
            n_molecules=request.n_molecules,
            pocket_center=request.pocket_center,
            pocket_radius=request.pocket_radius,
            task=request.task,
        )

        # Run PocketXMol via subprocess
        _run_sample_use(config_path, outdir)

        # Read SDF outputs → SMILES
        molecules = _parse_sdf_outputs(outdir)

        gen_time = time.time() - start_time
        return GenerateResponse(
            status="success",
            molecules=molecules,
            generation_time_seconds=round(gen_time, 1),
        )

    except Exception as e:
        gen_time = time.time() - start_time
        return GenerateResponse(
            status="error",
            molecules=[],
            generation_time_seconds=round(gen_time, 1),
            message=str(e),
        )

    finally:
        # Clean up PDB input
        try:
            os.unlink(pdb_path)
        except Exception:
            pass
        # Clean up output dir (optional — comment out to keep for debugging)
        # if outdir:
        #     import shutil
        #     shutil.rmtree(outdir, ignore_errors=True)


def _write_task_config(
    pdb_path: str,
    outdir: str,
    n_molecules: int,
    pocket_center: list[float],
    pocket_radius: float,
    task: str,
) -> str:
    """Generate a PocketXMol task YAML config dynamically.

    For de novo molecule design (no ligand input), use sbdd_simple template.
    For docking (with known ligand), use dock_smallmol template.
    """
    # Use sbdd_simple for de novo generation (no input ligand required)
    template_name = "sbdd_simple" if task in ("dock_smallmol", "sbdd") else task
    pxm_config_dir = os.path.join(_PXM_DIR, "configs", "sample", "examples")
    template_path = os.path.join(pxm_config_dir, f"{template_name}.yml")

    # Load template if exists, otherwise use minimal config
    if os.path.exists(template_path):
        with open(template_path) as f:
            config = yaml.safe_load(f)
    else:
        config = {}

    # Override with request parameters
    config.setdefault("sample", {})
    config["sample"]["num_mols"] = n_molecules
    config["sample"]["seed"] = 2024

    config.setdefault("data", {})
    config["data"]["protein_path"] = pdb_path
    config["data"]["is_pep"] = False
    # Remove input_ligand — we're doing de novo generation
    config["data"].pop("input_ligand", None)

    config["data"].setdefault("pocket_args", {})
    config["data"]["pocket_args"]["pocket_coord"] = pocket_center
    config["data"]["pocket_args"]["radius"] = pocket_radius

    config["data"].setdefault("pocmol_args", {})
    config["data"]["pocmol_args"]["data_id"] = f"api_{int(time.time())}"

    config.setdefault("transforms", {})
    config["transforms"].setdefault("featurizer_pocket", {})
    config["transforms"]["featurizer_pocket"]["center"] = pocket_center

    # Ensure variable_mol_size for sbdd (controls generated molecule size)
    if "variable_mol_size" not in config["transforms"]:
        config["transforms"]["variable_mol_size"] = {
            "name": "variable_mol_size",
            "num_atoms_distri": {
                "strategy": "mol_atoms_based",
                "mean": {"coef": 0, "bias": 28},
                "std": {"coef": 0, "bias": 2},
                "min": 5,
            },
        }

    # Set task type and noise
    config.setdefault("task", {})
    config["task"]["name"] = "sbdd"
    config["task"]["transform"] = {"name": "sbdd"}

    config.setdefault("noise", {})
    config["noise"]["name"] = "sbdd"

    # Write to temp file
    config_out = os.path.join(outdir, "task_config.yml")
    with open(config_out, "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    return config_out


def _run_sample_use(config_path: str, outdir: str):
    """Run PocketXMol's sample_use.py as a subprocess."""
    sample_script = os.path.join(_PXM_DIR, "scripts", "sample_use.py")
    model_config = os.path.join(_PXM_DIR, "configs", "sample", "pxm.yml")

    cmd = [
        sys.executable, sample_script,
        "--config_task", config_path,
        "--config_model", model_config,
        "--outdir", outdir,
        "--device", _DEVICE,
    ]

    print(f"[PXM] Running: {' '.join(cmd)}", flush=True)
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=1800,  # 30 min max
        cwd=_PXM_DIR,  # Run from PocketXMol directory (imports depend on this)
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"PocketXMol failed (exit {result.returncode}):\n"
            f"STDOUT: {result.stdout[-1000:]}\n"
            f"STDERR: {result.stderr[-1000:]}"
        )

    print(f"[PXM] Completed successfully", flush=True)


def _parse_sdf_outputs(outdir: str) -> list[MoleculeResult]:
    """Read PocketXMol outputs: gen_info.csv for SMILES + main SDFs for properties."""
    from rdkit import Chem
    from rdkit.Chem import QED, Descriptors

    molecules = []
    out_path = Path(outdir)

    # Find the run subdirectory (named like "task_config_pxm_*")
    run_dirs = sorted(out_path.glob("task_config_pxm_*"))
    if not run_dirs:
        return molecules
    run_dir = run_dirs[0]

    # Strategy 1: Read gen_info.csv for SMILES (already reconstructed by PocketXMol)
    gen_csv = run_dir / "gen_info.csv"
    smiles_from_csv: dict[str, dict] = {}  # smiles -> row data
    if gen_csv.exists():
        import csv
        with open(gen_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                smi = row.get("smiles", "").strip()
                if smi and Chem.MolFromSmiles(smi) is not None:
                    smiles_from_csv[smi] = {
                        "cfd": float(row.get("cfd_traj", 0) or 0),
                        "cfd_pos": float(row.get("cfd_pos", 0) or 0),
                        "filename": row.get("filename", ""),
                    }

    # If we have gen_info.csv SMILES, use those (most reliable)
    if smiles_from_csv:
        for smi, meta in smiles_from_csv.items():
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue

            qed_val = mw_val = logp_val = sa_val = None
            try:
                qed_val = round(QED.qed(mol), 3)
                mw_val = round(Descriptors.MolWt(mol), 1)
                logp_val = round(Descriptors.MolLogP(mol), 2)
            except Exception:
                pass
            try:
                from rdkit.Chem import RDConfig
                sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
                if sa_path not in sys.path:
                    sys.path.append(sa_path)
                import sascorer
                sa_val = round(sascorer.calculateScore(mol), 2)
            except Exception:
                pass

            molecules.append(MoleculeResult(
                smiles=smi,
                score=round(meta["cfd"], 3),
                qed=qed_val,
                mw=mw_val,
                logp=logp_val,
                sa_score=sa_val,
            ))
        return molecules

    # Strategy 2: Fallback — parse main SDF files only (not raw intermediates)
    sdf_dir = run_dir / f"{run_dir.name}_SDF"
    if sdf_dir.exists():
        for sdf_file in sorted(sdf_dir.glob("[0-9]*.sdf")):
            supplier = Chem.SDMolSupplier(str(sdf_file))
            for mol in supplier:
                if mol is None:
                    continue
                smiles = Chem.MolToSmiles(mol, canonical=True)
                if not smiles:
                    continue

                qed_val = mw_val = logp_val = sa_val = None
                try:
                    qed_val = round(QED.qed(mol), 3)
                    mw_val = round(Descriptors.MolWt(mol), 1)
                    logp_val = round(Descriptors.MolLogP(mol), 2)
                except Exception:
                    pass
                try:
                    from rdkit.Chem import RDConfig
                    sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
                    if sa_path not in sys.path:
                        sys.path.append(sa_path)
                    import sascorer
                    sa_val = round(sascorer.calculateScore(mol), 2)
                except Exception:
                    pass

                score = 0.0
                for prop_name in ("confidence", "score", "cfd_traj"):
                    if mol.HasProp(prop_name):
                        try:
                            score = float(mol.GetProp(prop_name))
                            break
                        except Exception:
                            pass

                molecules.append(MoleculeResult(
                    smiles=smiles,
                    score=score,
                    qed=qed_val,
                    mw=mw_val,
                    logp=logp_val,
                    sa_score=sa_val,
                ))

    return molecules


def main():
    global _PXM_DIR, _DEVICE, _CHECKPOINT_EXISTS

    parser = argparse.ArgumentParser(description="PocketXMol GPU Inference Server")
    parser.add_argument("--pxm-dir", required=True,
                        help="Path to PocketXMol installation (e.g., /root/PocketXMol)")
    parser.add_argument("--host", default="0.0.0.0", help="Server host (default: 0.0.0.0)")
    parser.add_argument("--port", type=int, default=8000, help="Server port (default: 8000)")
    parser.add_argument("--device", default="cuda:0", help="Device for inference (default: cuda:0)")
    args = parser.parse_args()

    _PXM_DIR = os.path.abspath(args.pxm_dir)
    _DEVICE = args.device

    # Verify PocketXMol installation
    sample_script = os.path.join(_PXM_DIR, "scripts", "sample_use.py")
    if not os.path.exists(sample_script):
        print(f"[ERROR] sample_use.py not found at {sample_script}")
        sys.exit(1)

    # Check model checkpoint (pxm_use for newer releases, pxm for older)
    ckpt_path = os.path.join(
        _PXM_DIR, "data", "trained_models", "pxm_use", "checkpoints", "pocketxmol.ckpt"
    )
    if not os.path.exists(ckpt_path):
        ckpt_path = os.path.join(
            _PXM_DIR, "data", "trained_models", "pxm", "checkpoints", "pocketxmol.ckpt"
        )
    _CHECKPOINT_EXISTS = os.path.exists(ckpt_path)
    if _CHECKPOINT_EXISTS:
        print(f"[OK] Checkpoint found: {ckpt_path}")
    else:
        print(f"[WARN] Checkpoint not found: {ckpt_path}")
        print(f"       Download from https://zenodo.org/records/17801271")
        print(f"       Extract: tar -zxvf model_weights.tar.gz -C {_PXM_DIR}")
        print(f"       Server will start but /generate will return 503 until checkpoint is available")

    import uvicorn
    print(f"[OK] Starting server on {args.host}:{args.port}")
    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
