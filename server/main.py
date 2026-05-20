#!/usr/bin/env python3
"""PocketXMol GPU Inference Server.

FastAPI service that exposes PocketXMol molecule generation via HTTP API.
Deploy on a GPU server and set DIFFUSION_API_URL on the Agent machine.

Usage:
    pip install -r requirements.txt
    python main.py --host 0.0.0.0 --port 8000 --model-dir /path/to/model_weights

Endpoints:
    GET  /health  - Health check
    POST /generate - Generate molecules for a protein pocket
"""
import argparse
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

app = FastAPI(title="PocketXMol Inference Server")

# Global model reference (loaded once at startup)
_model = None
_device = None
_model_dir = None


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


@app.get("/health", response_model=HealthResponse)
def health_check():
    import torch
    return HealthResponse(
        status="ok" if _model is not None else "model_not_loaded",
        model_loaded=_model is not None,
        gpu_available=torch.cuda.is_available(),
        device=_device or "",
    )


@app.post("/generate", response_model=GenerateResponse)
def generate_molecules(request: GenerateRequest):
    if _model is None:
        raise HTTPException(status_code=503, detail="Model not loaded")

    start_time = time.time()

    try:
        # Write PDB content to temp file (PocketXMol reads from file)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False, prefix="pxm_input_"
        ) as f:
            f.write(request.pdb_content)
            pdb_path = f.name

        try:
            molecules = _run_pocketxmol(
                pdb_path=pdb_path,
                n_molecules=request.n_molecules,
                pocket_center=request.pocket_center,
                pocket_radius=request.pocket_radius,
                task=request.task,
            )
        finally:
            os.unlink(pdb_path)

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


def _run_pocketxmol(
    pdb_path: str,
    n_molecules: int,
    pocket_center: list[float],
    pocket_radius: float,
    task: str,
) -> list[MoleculeResult]:
    """Run PocketXMol inference and convert results to SMILES."""
    import torch
    from rdkit import Chem
    from rdkit.Chem import QED, Descriptors

    # PocketXMol inference
    with tempfile.TemporaryDirectory(prefix="pxm_output_") as outdir:
        _sample_molecules(
            model=_model,
            pdb_path=pdb_path,
            outdir=outdir,
            n_molecules=n_molecules,
            pocket_center=pocket_center,
            pocket_radius=pocket_radius,
            device=_device,
        )

        # Read generated SDF files and convert to SMILES
        molecules = []
        sdf_dir = Path(outdir)
        for sdf_file in sorted(sdf_dir.glob("**/*.sdf")):
            supplier = Chem.SDMolSupplier(str(sdf_file))
            for mol in supplier:
                if mol is None:
                    continue
                smiles = Chem.MolToSmiles(mol, canonical=True)
                if not smiles:
                    continue

                # Compute drug-likeness properties
                qed_val = None
                mw_val = None
                logp_val = None
                sa_val = None
                try:
                    qed_val = round(QED.qed(mol), 3)
                    mw_val = round(Descriptors.MolWt(mol), 1)
                    logp_val = round(Descriptors.MolLogP(mol), 2)
                except Exception:
                    pass
                try:
                    from rdkit.Chem import RDConfig
                    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
                    import sascorer
                    sa_val = round(sascorer.calculateScore(mol), 2)
                except Exception:
                    pass

                # Extract score from SDF properties if available
                score = 0.0
                if mol.HasProp("score"):
                    try:
                        score = float(mol.GetProp("score"))
                    except Exception:
                        pass
                if mol.HasProp("confidence"):
                    try:
                        score = float(mol.GetProp("confidence"))
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


def _sample_molecules(
    model,
    pdb_path: str,
    outdir: str,
    n_molecules: int,
    pocket_center: list[float],
    pocket_radius: float,
    device: str,
):
    """Run PocketXMol sampling.

    This function calls PocketXMol's inference API. The exact call pattern
    depends on PocketXMol's version — adjust imports and calls as needed.
    """
    import torch

    # Try the high-level sample_use.py API first
    try:
        from PocketXMol.scripts.sample_use import sample as pxm_sample
        pxm_sample(
            model=model,
            pdb_path=pdb_path,
            outdir=outdir,
            n_molecules=n_molecules,
            pocket_center=pocket_center,
            pocket_radius=pocket_radius,
            device=device,
        )
        return
    except (ImportError, AttributeError):
        pass

    # Fallback: call sample_use.py as a subprocess
    import subprocess
    config_path = _find_config_for_task("dock_smallmol")
    cmd = [
        sys.executable, "-m", "PocketXMol.scripts.sample_use",
        "--config_task", config_path,
        "--outdir", outdir,
        "--device", device,
        "--pdb_path", pdb_path,
        "--n_molecules", str(n_molecules),
        "--pocket_center", json.dumps(pocket_center),
        "--pocket_radius", str(pocket_radius),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if result.returncode != 0:
        raise RuntimeError(f"PocketXMol failed: {result.stderr[:500]}")


def _find_config_for_task(task: str) -> str:
    """Find the PocketXMol config file for a given task."""
    # Common config locations
    search_paths = [
        Path(_model_dir) / "configs" / "sample" / "examples" / f"{task}.yml",
        Path(_model_dir) / "configs" / f"{task}.yml",
        Path("configs/sample/examples") / f"{task}.yml",
    ]
    for p in search_paths:
        if p.exists():
            return str(p)
    # Return default path and let PocketXMol handle the error
    return f"configs/sample/examples/{task}.yml"


def load_model(model_dir: str, device: str):
    """Load PocketXMol model weights at startup."""
    global _model, _device, _model_dir

    _model_dir = model_dir
    _device = device

    try:
        import torch

        # Try loading via PocketXMol's API
        try:
            from PocketXMol.models import load_pretrained
            _model = load_pretrained(model_dir, device=device)
            print(f"[OK] Model loaded from {model_dir} on {device}")
            return
        except (ImportError, AttributeError):
            pass

        # Alternative: load checkpoint directly
        ckpt_path = Path(model_dir) / "model.pt"
        if not ckpt_path.exists():
            ckpt_path = Path(model_dir) / "checkpoint.pt"
        if not ckpt_path.exists():
            print(f"[WARN] No checkpoint found in {model_dir}, model will not be available")
            return

        _model = torch.load(str(ckpt_path), map_location=device)
        _model.eval()
        print(f"[OK] Checkpoint loaded from {ckpt_path} on {device}")

    except Exception as e:
        print(f"[ERROR] Failed to load model: {e}")
        _model = None


def main():
    parser = argparse.ArgumentParser(description="PocketXMol GPU Inference Server")
    parser.add_argument("--host", default="0.0.0.0", help="Server host (default: 0.0.0.0)")
    parser.add_argument("--port", type=int, default=8000, help="Server port (default: 8000)")
    parser.add_argument("--model-dir", required=True, help="Path to PocketXMol model weights directory")
    parser.add_argument("--device", default="cuda:0", help="Device for inference (default: cuda:0)")
    parser.add_argument("--reload", action="store_true", help="Enable auto-reload for development")
    args = parser.parse_args()

    load_model(args.model_dir, args.device)

    import uvicorn
    uvicorn.run(
        "main:app" if args.reload else app,
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


if __name__ == "__main__":
    main()
