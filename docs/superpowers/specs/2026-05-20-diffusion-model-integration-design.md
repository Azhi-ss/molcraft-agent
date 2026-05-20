# PocketXMol Diffusion Model Integration Design

## Context

MolCraft Agent currently uses RDKit scaffold mutation (`src/generator.py`) for molecule generation. This approach is blind to the protein pocket and explores chemical space randomly. PocketXMol is an atom-level generative foundation model that generates 3D molecules inside protein binding pockets, achieving SOTA on 11/13 benchmarks. Integrating it as a supplementary generator will provide pocket-aware initial seeds for the evolutionary pipeline.

## Architecture

```
Agent (Kimi SDK)
  │
  ├─ DiffusionGenerate 工具 ──HTTP POST──→ GPU Server
  │       ↓                                    │
  │  SMILES 列表 ←────JSON Response───────────┘
  │       ↓
  ├─ run_pipeline (双引擎模式)
  │    ├─ diffusion: 初始生成 → dock → evolve
  │    ├─ hybrid: 扩散+RDKit 混合
  │    └─ mutate: RDKit 变异 (现有默认)
  │
  └─ 降级: GPU 不可用时自动回退到 RDKit
```

## GPU Server: FastAPI Service

### Endpoints

| Endpoint | Method | Input | Output | Description |
|----------|--------|-------|--------|-------------|
| `/health` | GET | - | `{status, model_loaded, gpu_available}` | Health check |
| `/generate` | POST | PDB content + params | `[{smiles, score, qed, sa_score}...]` | Molecule generation |

### `/generate` Request Format

```json
{
  "pdb_content": "ATOM  ...",
  "n_molecules": 20,
  "pocket_center": [12.5, -3.2, 8.7],
  "pocket_radius": 10.0,
  "task": "dock_smallmol"
}
```

### `/generate` Response Format

```json
{
  "status": "success",
  "molecules": [
    {
      "smiles": "c1ccc(cc1)C(=O)N",
      "score": 0.85,
      "qed": 0.72,
      "sa_score": 3.1
    }
  ],
  "generation_time_seconds": 120.5
}
```

### Server Responsibilities

1. Receive PDB content (not file path — server is stateless between calls)
2. Write PDB to temp file, invoke PocketXMol
3. Convert output SDF → SMILES (server-side, using RDKit)
4. Compute QED and SA score (server-side, using RDKit)
5. Return JSON with SMILES list
6. Handle timeout: if inference exceeds `DIFFUSION_TIMEOUT`, return error

### Server Implementation

- FastAPI app in `server/` directory (not part of Agent's main codebase)
- PocketXMol loaded once at startup, reused for all requests
- Async inference with `asyncio.to_thread` to avoid blocking the event loop
- Config via environment variables: `MODEL_DIR`, `DEVICE`, `BATCH_SIZE`

## Agent-Side: New Tool

### `DiffusionGenerate` Tool

Registered in `agent.yaml` as `molcraft_agent.tools:DiffusionGenerate`.

```python
class DiffusionGenerateParams(BaseModel):
    pdb_path: str = Field(default="data/target.pdb")
    n_molecules: int = Field(default=20)
    pocket_center: list[float] | None = Field(default=None)

class DiffusionGenerate(CallableTool2):
    name = "diffusion_generate"
    description = (
        "用 PocketXMol 扩散模型生成口袋感知分子。需要 GPU 服务器。"
        "返回 SMILES 列表及药物性质，格式与 generate_molecules 兼容。"
    )
```

### Call Flow

1. Read PDB file content
2. POST to `DIFFUSION_API_URL/generate`
3. If timeout > 60s: poll `/health` every 30s until result ready
4. Return SMILES list in same format as `generate_molecules` output
5. If GPU server unreachable: return `ToolError` with fallback hint

### Tool Registration

Add to `agent.yaml`:

```yaml
- "molcraft_agent.tools:DiffusionGenerate"
```

## Pipeline Integration

### New `--generator` Parameter

In `tools/pipeline.py` and `RunPipeline` tool:

| Value | Behavior |
|-------|----------|
| `mutate` | Existing RDKit mutation (default, backward compatible) |
| `diffusion` | Pure PocketXMol generation |
| `hybrid` | PocketXMol seeds + RDKit diversity supplement |

### Hybrid Mode Detail

1. Generate `n_diffusion` molecules via PocketXMol (e.g., 20)
2. Generate `n_rdkit` molecules via RDKit (e.g., 30)
3. Merge into single candidate pool
4. Proceed with existing docking → evolution → synthesis flow
5. Rationale: diffusion provides high-quality pocket-aware seeds, RDKit provides diversity

### Fallback

- If `DIFFUSION_API_URL` is not set or server is unreachable:
  - `diffusion` mode: ToolError, Agent must switch to mutate
  - `hybrid` mode: Automatically degrades to pure RDKit
  - No agent crash or hang

## Configuration

In `src/config.py` or `.env`:

```python
DIFFUSION_API_URL = os.getenv("DIFFUSION_API_URL", "")  # e.g., http://192.168.1.100:8000
DIFFUSION_ENABLED = bool(DIFFUSION_API_URL)
DIFFUSION_TIMEOUT = int(os.getenv("DIFFUSION_TIMEOUT", "1800"))  # 30 min
DIFFUSION_DEFAULT_N = int(os.getenv("DIFFUSION_DEFAULT_N", "20"))
```

## Files to Create/Modify

| File | Action | Description |
|------|--------|-------------|
| `server/main.py` | Create | FastAPI service for PocketXMol |
| `server/requirements.txt` | Create | GPU server dependencies |
| `molcraft_agent/tools.py` | Modify | Add `DiffusionGenerate` class |
| `agent.yaml` | Modify | Register new tool |
| `src/config.py` | Modify | Add diffusion config vars |
| `tools/pipeline.py` | Modify | Add `--generator` parameter |
| `program.md` | Modify | Update Agent instructions for diffusion usage |

## Out of Scope

- Training or fine-tuning PocketXMol (use pretrained weights only)
- Replacing RDKit generator entirely (dual-engine approach)
- Real-time streaming of generation progress (polling is sufficient)
- Multiple GPU server load balancing (single server for now)
