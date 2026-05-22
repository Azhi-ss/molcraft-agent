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
_START_TIME: float = time.time()
_PXM_DIR: str = ""
_DEVICE: str = "cuda:0"
_CHECKPOINT_EXISTS: bool = False
_jobs: dict[str, dict] = {}
_jobs_lock = threading.Lock()


def _gpu_info() -> tuple[int, int, int]:
    """Return (memory_used_mb, memory_total_mb, utilization_pct)."""
    try:
        import subprocess
        r = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.used,memory.total,utilization.gpu",
             "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=5,
        )
        parts = r.stdout.strip().split(",")
        return int(parts[0].strip()), int(parts[1].strip()), int(parts[2].strip())
    except Exception:
        return 0, 0, 0

# ── Models ──


class GenerateRequest(BaseModel):
    pdb_content: str = Field(description="Protein PDB file content")
    n_molecules: int = Field(default=20, ge=1, le=200, description="Number of molecules to generate (1-200)")
    pocket_center: list[float] = Field(
        default=[0.0, 0.0, 0.0],
        description="Pocket center coordinates [x, y, z]",
    )
    pocket_radius: float = Field(default=15.0, description="Pocket radius in Angstroms")
    task: str = Field(default="sbdd", description="Task type: sbdd (de novo AR refine, recommended), sbdd_simple (fast)")


class OptimizeRequest(BaseModel):
    pdb_content: str = Field(description="Protein PDB file content")
    smiles: str = Field(description="Seed molecule SMILES to optimize")
    n_variants: int = Field(default=10, description="Number of optimized variants to generate")
    pocket_center: list[float] | None = Field(default=None, description="Pocket center, derived from SMILES if not set")


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
    binding_energy: Optional[float] = None


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
    gpu_memory_used_mb: int = 0
    gpu_memory_total_mb: int = 0
    gpu_utilization_pct: int = 0
    device: str = ""
    pxm_dir: str = ""


class StatsResponse(BaseModel):
    uptime_seconds: float
    total_jobs_completed: int
    total_jobs_failed: int
    total_molecules_generated: int
    gpu_memory_used_mb: int
    gpu_memory_total_mb: int
    gpu_utilization_pct: int


# ── Endpoints ──


@app.get("/health", response_model=HealthResponse)
def health_check():
    import torch
    with _jobs_lock:
        active = sum(1 for j in _jobs.values() if j["status"] == "running")
    gpu_mem_used, gpu_mem_total, gpu_util = _gpu_info()
    return HealthResponse(
        status="ok" if _CHECKPOINT_EXISTS else "checkpoint_missing",
        model_loaded=_CHECKPOINT_EXISTS,
        gpu_available=torch.cuda.is_available(),
        active_jobs=active,
        gpu_memory_used_mb=gpu_mem_used,
        gpu_memory_total_mb=gpu_mem_total,
        gpu_utilization_pct=gpu_util,
        device=_DEVICE,
        pxm_dir=_PXM_DIR,
    )


@app.get("/stats", response_model=StatsResponse)
def get_stats():
    with _jobs_lock:
        completed = sum(1 for j in _jobs.values() if j["status"] == "completed")
        failed = sum(1 for j in _jobs.values() if j["status"] == "failed")
        total_mols = sum(len(j.get("result") or []) for j in _jobs.values() if j["status"] == "completed")

    gpu_mem_used, gpu_mem_total, gpu_util = _gpu_info()
    return StatsResponse(
        uptime_seconds=round(time.time() - _START_TIME, 1),
        total_jobs_completed=completed,
        total_jobs_failed=failed,
        total_molecules_generated=total_mols,
        gpu_memory_used_mb=gpu_mem_used,
        gpu_memory_total_mb=gpu_mem_total,
        gpu_utilization_pct=gpu_util,
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


@app.post("/optimize", response_model=GenerateAccepted, status_code=202)
def submit_optimization(request: OptimizeRequest):
    """Optimize a seed molecule around the pocket (opt_mol task)."""
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
        target=_run_opt_job,
        args=(job_id, request),
        daemon=True,
    )
    thread.start()

    with _jobs_lock:
        _jobs[job_id]["thread"] = thread

    print(f"[Job {job_id}] Optimize: smiles={request.smiles[:40]}... n={request.n_variants}", flush=True)
    return GenerateAccepted(
        status="running",
        job_id=job_id,
        message=f"Optimization accepted, {request.n_variants} variants requested",
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
        if outdir:
            import shutil
            shutil.rmtree(outdir, ignore_errors=True)


def _run_opt_job(job_id: str, request: OptimizeRequest):
    """Run PocketXMol opt_mol optimization in background."""
    start = time.time()
    pdb_path = None
    ligand_path = None
    outdir = None

    try:
        # Write PDB
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False, prefix="pxm_input_"
        ) as f:
            f.write(request.pdb_content)
            pdb_path = f.name

        # Write ligand SMILES → 3D SDF
        from rdkit import Chem
        ligand_path = tempfile.mktemp(suffix=".sdf", prefix="pxm_ligand_")
        mol = Chem.MolFromSmiles(request.smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {request.smiles}")
        mol = Chem.AddHs(mol)
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(mol, randomSeed=2024)
        AllChem.MMFFOptimizeMolecule(mol)
        writer = Chem.SDWriter(ligand_path)
        writer.write(mol)
        writer.close()

        # Pocket center from ligand if not set
        if request.pocket_center:
            pocket_center = request.pocket_center
        else:
            conf = mol.GetConformer()
            coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            import numpy as np
            center = np.mean([[a.x, a.y, a.z] for a in coords], axis=0)
            pocket_center = [round(float(c), 2) for c in center]

        outdir = tempfile.mkdtemp(prefix="pxm_output_")

        config_path = _write_opt_config(
            pdb_path=pdb_path,
            ligand_path=ligand_path,
            outdir=outdir,
            n_variants=request.n_variants,
            pocket_center=pocket_center,
        )

        _run_sample_use(config_path, outdir)
        molecules = _parse_sdf_outputs(outdir)

        elapsed = time.time() - start
        with _jobs_lock:
            _jobs[job_id]["status"] = "completed"
            _jobs[job_id]["result"] = molecules
            _jobs[job_id]["finished_at"] = time.time()

        print(f"[Job {job_id}] Opt done: {len(molecules)} variants in {elapsed:.0f}s", flush=True)

    except Exception as e:
        elapsed = time.time() - start
        with _jobs_lock:
            _jobs[job_id]["status"] = "failed"
            _jobs[job_id]["error"] = str(e)
        print(f"[Job {job_id}] Opt failed ({elapsed:.0f}s): {e}", flush=True)

    finally:
        for path in [pdb_path, ligand_path]:
            if path:
                try:
                    os.unlink(path)
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
    # sbdd_simple is the working de novo config
    # sbdd (AR refine) needs train prior files not available — use only if available
    if task == "sbdd":
        template_name = "sbdd_simple"  # Fallback to working config
    elif task in ("dock_smallmol", "sbdd_simple"):
        template_name = "sbdd_simple"
    else:
        template_name = task
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


def _write_opt_config(
    pdb_path: str,
    ligand_path: str,
    outdir: str,
    n_variants: int,
    pocket_center: list[float],
) -> str:
    """Generate opt_mol task config — optimizes input ligand around pocket."""
    pxm_config_dir = os.path.join(_PXM_DIR, "configs", "sample", "examples")
    template_path = os.path.join(pxm_config_dir, "opt_mol.yml")

    if os.path.exists(template_path):
        with open(template_path) as f:
            config = yaml.safe_load(f)
    else:
        config = {}

    config.setdefault("sample", {})
    config["sample"]["num_mols"] = n_variants
    config["sample"]["seed"] = 2024

    config.setdefault("data", {})
    config["data"]["protein_path"] = pdb_path
    config["data"]["input_ligand"] = ligand_path
    config["data"]["is_pep"] = False
    config["data"].setdefault("pocket_args", {})
    config["data"]["pocket_args"]["radius"] = 10

    config.setdefault("transforms", {})
    config["transforms"].setdefault("featurizer", {})
    config["transforms"]["featurizer"]["mol_as_pocket_center"] = True

    # Ensure variable_mol_size (slightly larger for optimization)
    if "variable_mol_size" not in config["transforms"]:
        config["transforms"]["variable_mol_size"] = {
            "name": "variable_mol_size",
            "num_atoms_distri": {
                "strategy": "mol_atoms_based",
                "mean": {"coef": 0, "bias": 38},
                "std": {"coef": 0, "bias": 3},
                "min": 5,
            },
        }

    config.setdefault("task", {})
    config["task"]["name"] = "sbdd"
    config["task"]["transform"] = {"name": "sbdd"}

    config.setdefault("noise", {})
    config["noise"]["name"] = "sbdd"
    config["noise"]["num_steps"] = 50
    config["noise"]["init_step"] = 0.5  # Start closer to input → more similar output

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


# ── Generate + Dock endpoint (diffusion → docking in one call) ──

_RECEPTOR_PDBQT: str | None = None


def _prepare_receptor(pdb_content: str) -> str:
    """Prepare receptor PDBQT from PDB content. Cached after first call. Thread-safe."""
    global _RECEPTOR_PDBQT
    with _jobs_lock:
        if _RECEPTOR_PDBQT and os.path.exists(_RECEPTOR_PDBQT):
            return _RECEPTOR_PDBQT

    # Build receptor outside lock (only one thread reaches here)
    pdb_path = tempfile.mktemp(suffix=".pdb", prefix="receptor_")
    with open(pdb_path, "w") as f:
        f.write(pdb_content)

    pdbqt_path = tempfile.mktemp(suffix=".pdbqt", prefix="receptor_")
    subprocess.run(
        ["obabel", pdb_path, "-O", pdbqt_path, "-xr"],
        capture_output=True, text=True, timeout=60,
        check=True,
    )
    os.unlink(pdb_path)

    with _jobs_lock:
        if _RECEPTOR_PDBQT is None:
            _RECEPTOR_PDBQT = pdbqt_path
        else:
            os.unlink(pdbqt_path)  # Another thread finished first, discard ours
            pdbqt_path = _RECEPTOR_PDBQT
    return pdbqt_path


def _dock_single_from_sdf(sdf_file: str, receptor_pdbqt: str,
                          center: list[float], size: list[float]) -> float | None:
    """Dock a molecule from SDF via obabel → PDBQT → Vina. Returns binding_energy."""
    pdbqt_path = tempfile.mktemp(suffix=".pdbqt", prefix="ligand_")

    try:
        # obabel: SDF → PDBQT (preserves PocketXMol's 3D pocket coordinates)
        subprocess.run(
            ["obabel", sdf_file, "-O", pdbqt_path],
            capture_output=True, text=True, timeout=30, check=True,
        )
    except subprocess.CalledProcessError:
        return None

    try:
        out_path = tempfile.mktemp(suffix=".pdbqt", prefix="vina_out_")
        result = subprocess.run(
            [
                "vina",
                "--receptor", receptor_pdbqt,
                "--ligand", pdbqt_path,
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(size[0]),
                "--size_y", str(size[1]),
                "--size_z", str(size[2]),
                "--out", out_path,
                "--exhaustiveness", "8",
            ],
            capture_output=True, text=True, timeout=120,
        )
    except subprocess.TimeoutExpired:
        return None
    finally:
        os.unlink(pdbqt_path)
        if os.path.exists(out_path):
            os.unlink(out_path)

    for line in result.stdout.split("\n"):
        if line.strip().startswith("1 "):
            try:
                return round(float(line.strip().split()[1]), 3)
            except (IndexError, ValueError):
                pass
    return None


@app.post("/generate-and-dock", response_model=GenerateAccepted, status_code=202)
def submit_generate_and_dock(request: GenerateRequest):
    """Generate molecules via PocketXMol + dock on GPU server. Returns with binding_energy."""
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
        target=_run_gen_and_dock_job,
        args=(job_id, request),
        daemon=True,
    )
    thread.start()

    with _jobs_lock:
        _jobs[job_id]["thread"] = thread

    print(f"[Job {job_id}] Gen+Dock: n={request.n_molecules}", flush=True)
    return GenerateAccepted(
        status="running",
        job_id=job_id,
        message=f"Gen+ dock accepted, {request.n_molecules} molecules requested",
    )


def _run_gen_and_dock_job(job_id: str, request: GenerateRequest):
    """Generate → dock pipeline on GPU server."""
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

        # Stage 1: Generate
        config_path = _write_task_config(
            pdb_path=pdb_path, outdir=outdir,
            n_molecules=request.n_molecules,
            pocket_center=request.pocket_center,
            pocket_radius=request.pocket_radius,
            task=request.task,
        )
        _run_sample_use(config_path, outdir)

        # Find main SDF directory
        run_dirs = sorted(Path(outdir).glob("task_config_pxm_*"))
        if not run_dirs:
            raise RuntimeError("No PocketXMol output directory found")
        sdf_dir = run_dirs[0] / f"{run_dirs[0].name}_SDF"

        # Stage 2: Dock main SDFs directly (preserves PocketXMol 3D coords)
        receptor_pdbqt = _prepare_receptor(request.pdb_content)
        box_size = [25.0, 25.0, 25.0]

        be_by_filename: dict[str, float] = {}
        if sdf_dir.exists():
            for sdf_file in sorted(sdf_dir.glob("[0-9]*.sdf")):
                if "bad" in sdf_file.name:
                    continue
                be = _dock_single_from_sdf(
                    str(sdf_file), receptor_pdbqt,
                    request.pocket_center, box_size,
                )
                if be is not None:
                    be_by_filename[sdf_file.name] = be
                print(f"  [{sdf_file.name}] BE={be}", flush=True)

        # Stage 3: Build results from gen_info.csv + docking scores
        from rdkit import Chem
        from rdkit.Chem import QED, Descriptors

        molecules: list[MoleculeResult] = []
        gen_csv = run_dirs[0] / "gen_info.csv"
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

                    fname = row.get("filename", "")
                    be = be_by_filename.get(fname)

                    molecules.append(MoleculeResult(
                        smiles=smi,
                        score=round(cfd, 3),
                        qed=qed_val,
                        mw=mw_val,
                        logp=logp_val,
                        sa_score=sa_val,
                        binding_energy=be,
                    ))

        docked_count = sum(1 for m in molecules if m.binding_energy is not None)
        elapsed = time.time() - start
        with _jobs_lock:
            _jobs[job_id]["status"] = "completed"
            _jobs[job_id]["result"] = [m.model_dump() for m in molecules]
            _jobs[job_id]["finished_at"] = time.time()

        print(f"[Job {job_id}] Gen+Dock done: {len(molecules)} mols, "
              f"{docked_count} docked in {elapsed:.0f}s", flush=True)

    except Exception as e:
        elapsed = time.time() - start
        with _jobs_lock:
            _jobs[job_id]["status"] = "failed"
            _jobs[job_id]["error"] = str(e)
        print(f"[Job {job_id}] Gen+Dock failed ({elapsed:.0f}s): {e}", flush=True)

    finally:
        if pdb_path:
            try:
                os.unlink(pdb_path)
            except Exception:
                pass
        if outdir:
            import shutil
            shutil.rmtree(outdir, ignore_errors=True)


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
