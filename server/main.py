#!/usr/bin/env python3
"""PocketXMol GPU Inference Server — Async Job Mode.

FastAPI service with async job queue. POST /generate creates a job and returns
immediately; Agent polls GET /job/{id} for status and GET /job/{id}/result for output.

Usage:
    cd /root/PocketXMol
    pip install fastapi uvicorn pyyaml
    python server/main.py --pxm-dir /root/PocketXMol --port 8000

Endpoints:
    GET  /health          - Health + GPU check
    POST /generate        - Submit job → {job_id, status: "running"}
    GET  /job/{job_id}    - Job status {status, elapsed_s, message}
    GET  /job/{job_id}/result - Final molecules (or 404/503 if not ready)
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
import uuid
import yaml
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

app = FastAPI(title="PocketXMol Inference Server")

# ── Global state ──
_PXM_DIR: str = ""
_DEVICE: str = "cuda:0"
_CHECKPOINT_EXISTS: bool = False
_jobs: dict[str, dict] = {}
_jobs_lock = threading.Lock()

# ── Models ──


class GenerateRequest(BaseModel):
    pdb_content: str = Field(description="Protein PDB file content")
    n_molecules: int = Field(default=20, description="Number of molecules to generate")
    pocket_center: list[float] = Field(
        default=[0.0, 0.0, 0.0],
        description="Pocket center coordinates [x, y, z]",
    )
    pocket_radius: float = Field(default=15.0, description="Pocket radius in Angstroms")
    task: str = Field(default="dock_smallmol", description="Task type")


class GenerateAccepted(BaseModel):
    status: str
    job_id: str
    message: str = ""


class MoleculeResult(BaseModel):
    smiles: str
    score: float = 0.0
    qed: Optional[float] = None
    mw: Optional[float] = None
    logp: Optional[float] = None
    sa_score: Optional[float] = None


class JobStatus(BaseModel):
    status: str  # "running" | "completed" | "failed"
    elapsed_seconds: float
    molecule_count: int = 0
    message: str = ""


class JobResult(BaseModel):
    status: str
    molecules: list[MoleculeResult]
    generation_time_seconds: float


class HealthResponse(BaseModel):
    status: str
    model_loaded: bool
    gpu_available: bool
    active_jobs: int
    device: str = ""
    pxm_dir: str = ""


# ── Endpoints ──


@app.get("/health", response_model=HealthResponse)
def health_check():
    import torch
    with _jobs_lock:
        active = sum(1 for j in _jobs.values() if j["status"] == "running")
    return HealthResponse(
        status="ok" if _CHECKPOINT_EXISTS else "checkpoint_missing",
        model_loaded=_CHECKPOINT_EXISTS,
        gpu_available=torch.cuda.is_available(),
        active_jobs=active,
        device=_DEVICE,
        pxm_dir=_PXM_DIR,
    )


@app.post("/generate", response_model=GenerateAccepted, status_code=202)
def submit_generation(request: GenerateRequest):
    if not _CHECKPOINT_EXISTS:
        raise HTTPException(status_code=503, detail="Model checkpoint not found")

    job_id = uuid.uuid4().hex[:12]

    with _jobs_lock:
        _jobs[job_id] = {
            "status": "running",
            "created_at": time.time(),
            "result": None,
            "error": None,
            "thread": None,
        }

    thread = threading.Thread(
        target=_run_job,
        args=(job_id, request),
        daemon=True,
    )
    thread.start()

    with _jobs_lock:
        _jobs[job_id]["thread"] = thread

    print(f"[Job {job_id}] Accepted: n={request.n_molecules}", flush=True)
    return GenerateAccepted(
        status="running",
        job_id=job_id,
        message=f"Job accepted, {request.n_molecules} molecules requested",
    )


@app.get("/job/{job_id}", response_model=JobStatus)
def get_job_status(job_id: str):
    with _jobs_lock:
        job = _jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    elapsed = time.time() - job["created_at"]
    mol_count = len(job["result"]) if job["result"] else 0
    return JobStatus(
        status=job["status"],
        elapsed_seconds=round(elapsed, 1),
        molecule_count=mol_count,
        message=job["error"] or "",
    )


@app.get("/job/{job_id}/result", response_model=JobResult)
def get_job_result(job_id: str):
    with _jobs_lock:
        job = _jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    if job["status"] == "running":
        elapsed = time.time() - job["created_at"]
        raise HTTPException(
            status_code=202,
            detail=f"Job still running ({elapsed:.0f}s elapsed)",
        )
    if job["status"] == "failed":
        raise HTTPException(status_code=500, detail=job["error"] or "Unknown error")

    elapsed = job.get("finished_at", time.time()) - job["created_at"]
    return JobResult(
        status="completed",
        molecules=job["result"] or [],
        generation_time_seconds=round(elapsed, 1),
    )


# ── Background job runner ──


def _run_job(job_id: str, request: GenerateRequest):
    """Run PocketXMol inference in background, update job on completion."""
    start = time.time()
    pdb_path = None
    outdir = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False, prefix="pxm_input_"
        ) as f:
            f.write(request.pdb_content)
            pdb_path = f.name

        outdir = tempfile.mkdtemp(prefix="pxm_output_")

        config_path = _write_task_config(
            pdb_path=pdb_path,
            outdir=outdir,
            n_molecules=request.n_molecules,
            pocket_center=request.pocket_center,
            pocket_radius=request.pocket_radius,
            task=request.task,
        )

        _run_sample_use(config_path, outdir)
        molecules = _parse_sdf_outputs(outdir)

        elapsed = time.time() - start
        with _jobs_lock:
            _jobs[job_id]["status"] = "completed"
            _jobs[job_id]["result"] = molecules
            _jobs[job_id]["finished_at"] = time.time()

        print(f"[Job {job_id}] Done: {len(molecules)} mols in {elapsed:.0f}s", flush=True)

    except Exception as e:
        elapsed = time.time() - start
        with _jobs_lock:
            _jobs[job_id]["status"] = "failed"
            _jobs[job_id]["error"] = str(e)

        print(f"[Job {job_id}] Failed ({elapsed:.0f}s): {e}", flush=True)

    finally:
        if pdb_path:
            try:
                os.unlink(pdb_path)
            except Exception:
                pass


# ── Background cleanup ──


def _cleanup_old_jobs():
    """Remove jobs older than 1 hour (daemon thread)."""
    while True:
        time.sleep(300)  # Every 5 minutes
        cutoff = time.time() - 3600
        with _jobs_lock:
            stale = [
                jid for jid, j in _jobs.items()
                if j["created_at"] < cutoff and j["status"] != "running"
            ]
            for jid in stale:
                del _jobs[jid]
        if stale:
            print(f"[Cleanup] Removed {len(stale)} old jobs", flush=True)


# ── PocketXMol integration (unchanged logic) ──


def _write_task_config(
    pdb_path: str,
    outdir: str,
    n_molecules: int,
    pocket_center: list[float],
    pocket_radius: float,
    task: str,
) -> str:
    template_name = "sbdd_simple" if task in ("dock_smallmol", "sbdd") else task
    pxm_config_dir = os.path.join(_PXM_DIR, "configs", "sample", "examples")
    template_path = os.path.join(pxm_config_dir, f"{template_name}.yml")

    if os.path.exists(template_path):
        with open(template_path) as f:
            config = yaml.safe_load(f)
    else:
        config = {}

    config.setdefault("sample", {})
    config["sample"]["num_mols"] = n_molecules
    config["sample"]["seed"] = 2024

    config.setdefault("data", {})
    config["data"]["protein_path"] = pdb_path
    config["data"]["is_pep"] = False
    config["data"].pop("input_ligand", None)

    config["data"].setdefault("pocket_args", {})
    config["data"]["pocket_args"]["pocket_coord"] = pocket_center
    config["data"]["pocket_args"]["radius"] = pocket_radius

    config["data"].setdefault("pocmol_args", {})
    config["data"]["pocmol_args"]["data_id"] = f"api_{int(time.time())}"

    config.setdefault("transforms", {})
    config["transforms"].setdefault("featurizer_pocket", {})
    config["transforms"]["featurizer_pocket"]["center"] = pocket_center

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

    config.setdefault("task", {})
    config["task"]["name"] = "sbdd"
    config["task"]["transform"] = {"name": "sbdd"}

    config.setdefault("noise", {})
    config["noise"]["name"] = "sbdd"

    config_out = os.path.join(outdir, "task_config.yml")
    with open(config_out, "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    return config_out


def _run_sample_use(config_path: str, outdir: str):
    sample_script = os.path.join(_PXM_DIR, "scripts", "sample_use.py")
    model_config = os.path.join(_PXM_DIR, "configs", "sample", "pxm.yml")

    cmd = [
        sys.executable, sample_script,
        "--config_task", config_path,
        "--config_model", model_config,
        "--outdir", outdir,
        "--device", _DEVICE,
    ]

    print(f"[PXM] Job start: {' '.join(cmd)}", flush=True)
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=1800,
        cwd=_PXM_DIR,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"PocketXMol failed (exit {result.returncode}):\n"
            f"STDOUT: {result.stdout[-1000:]}\n"
            f"STDERR: {result.stderr[-1000:]}"
        )

    print(f"[PXM] Completed", flush=True)


def _parse_sdf_outputs(outdir: str) -> list[MoleculeResult]:
    from rdkit import Chem
    from rdkit.Chem import QED, Descriptors

    molecules = []
    out_path = Path(outdir)

    run_dirs = sorted(out_path.glob("task_config_pxm_*"))
    if not run_dirs:
        return molecules
    run_dir = run_dirs[0]

    gen_csv = run_dir / "gen_info.csv"
    if gen_csv.exists():
        with open(gen_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                smi = row.get("smiles", "").strip()
                if not smi:
                    continue
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

                cfd = 0.0
                try:
                    cfd = float(row.get("cfd_traj", 0) or 0)
                except Exception:
                    pass

                molecules.append(MoleculeResult(
                    smiles=smi,
                    score=round(cfd, 3),
                    qed=qed_val,
                    mw=mw_val,
                    logp=logp_val,
                    sa_score=sa_val,
                ))
        return molecules

    # Fallback: parse SDF files
    sdf_dir = run_dir / f"{run_dir.name}_SDF"
    if sdf_dir.exists():
        for sdf_file in sorted(sdf_dir.glob("[0-9]*.sdf")):
            if "bad" in sdf_file.name:
                continue
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


# ── Entry point ──


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

    sample_script = os.path.join(_PXM_DIR, "scripts", "sample_use.py")
    if not os.path.exists(sample_script):
        print(f"[ERROR] sample_use.py not found at {sample_script}")
        sys.exit(1)

    ckpt_path = os.path.join(
        _PXM_DIR, "data", "trained_models", "pxm_use", "checkpoints", "pocketxmol.ckpt"
    )
    if not os.path.exists(ckpt_path):
        ckpt_path = os.path.join(
            _PXM_DIR, "data", "trained_models", "pxm", "checkpoints", "pocketxmol.ckpt"
        )
    _CHECKPOINT_EXISTS = os.path.exists(ckpt_path)
    if _CHECKPOINT_EXISTS:
        print(f"[OK] Checkpoint: {ckpt_path}")
    else:
        print(f"[WARN] Checkpoint not found: {ckpt_path}")

    # Start cleanup daemon
    threading.Thread(target=_cleanup_old_jobs, daemon=True).start()

    import uvicorn
    print(f"[OK] Starting on {args.host}:{args.port} (async jobs, device={_DEVICE})")
    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
