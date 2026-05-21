---
name: lkm-gaia-integration
description: Integrate Bohrium LKM search and Gaia formalize into molcraft-agent
---

# LKM + Gaia Integration Design

## Goal

Give molcraft-agent access to the Bohrium LKM scientific knowledge graph for literature search and reasoning-chain tracing, plus the ability to formalize findings into Gaia knowledge packages as iteration evidence.

## Approach

Tool + Skill hybrid (方案 A):

- **LKM search**: Python tool (`SearchLKM` in `tools.py`) for the 3 high-frequency endpoints. Fast, type-safe, consistent with existing tools.
- **Gaia formalize**: Kimi skills (SKILL.md) copied from `gaia-lkm-skills/`. Low-frequency, workflow-driven, flexible.
- **SearchWeb preserved**: Not replaced; LKM is priority for scientific claims, SearchWeb is fallback.

## Part 1: SearchLKM Tool

### Location

`molcraft_agent/tools.py` — new class `SearchLKM(CallableTool2)`

### Parameters

| Field | Type | Description |
|-------|------|-------------|
| `verb` | `str` | Sub-command: `search`, `reasoning`, or `papers_graph` |
| `query` | `str \| None` | Search query (for `search` verb) |
| `claim_id` | `str \| None` | Claim GCN ID (for `reasoning` verb) |
| `doi` | `str \| None` | Paper DOI (for `papers_graph` verb) |
| `title` | `str \| None` | Paper title (for `papers_graph` verb) |
| `top_k` | `int` | Max results for search (default 20) |
| `reasoning_only` | `bool` | Only return claims with reasoning chains (default false) |

### Implementation

- Internal HTTP calls using `urllib.request` (same as `lkm_search.py` logic, no Shell dependency)
- `accessKey` from `os.environ["LKM_ACCESS_KEY"]`; return `ToolError` with setup hint if missing
- Base URL: `https://open.bohrium.com/openapi/v1/lkm`
- Return full JSON response (preserve claim ids, source packages, reasoning chains for downstream formalize)
- Auto-log to `experiments.jsonl` via `append_experiment`

### Endpoints covered

| Verb | API Endpoint | Purpose |
|------|-------------|---------|
| `search` | `POST /search` | Search scientific claims by topic |
| `reasoning` | `GET /claims/{id}/reasoning` | Trace derivation chain of a claim |
| `papers_graph` | `POST /papers/graph` | Retrieve paper's full knowledge graph |

### Endpoints NOT covered (use Shell + lkm_search.py)

- `POST /reasoning/search` — low frequency
- `POST /variables/batch` — low frequency

## Part 2: Gaia Formalize Skills

### Skill directories to copy

From `/home/dministrator/Lab/clones/gaia-lkm-skills/skills/` → `.kimi/skills/lkm/`:

```
.kimi/skills/lkm/
├── orchestrator/SKILL.md
├── lkm-search/SKILL.md
├── lkm-search/scripts/lkm_search.py
├── lkm-search/references/api-contract.md
├── lkm-explorer/SKILL.md
├── lkm-explorer/references/
└── formalize/SKILL.md
```

### Skills NOT copied

- `evidence-subgraph` — irrelevant to molcraft
- `scholarly-synthesis` — irrelevant to molcraft
- `lkm-search-internal` — requires whitelist access

### Agent workflow

1. Use `SearchLKM` tool to search claims
2. Read `.kimi/skills/lkm/orchestrator/SKILL.md` for routing
3. Follow SOP via `lkm-explorer` or `formalize`
4. Use Shell to call `gaia` CLI for compile/check/infer

### Prerequisite

`gaia` CLI must be installed (`pip install gaia` or from SiliconEinstein/Gaia).

## Part 3: Config and Guidance Changes

### .env

Add:
```
LKM_ACCESS_KEY="4fc2fbe980fe452bbb0c8d09d80f056a"
```

### agent.yaml

Register new tool:
```yaml
- "molcraft_agent.tools:SearchLKM"
```

### program.md

In "阶段一：文献解析" section, add search priority:
```
文献检索优先级：
1. SearchLKM — 搜索 Bohrium LKM 科学知识图谱，获取结构化 claims 和推理链
2. SearchWeb — LKM 无结果时兜底，搜通用文献和博客
3. 当需要将 LKM 证据图谱化时，读取 .kimi/skills/lkm/orchestrator/SKILL.md 按 SOP 走 formalize
```

In "阶段二：瓶颈诊断" section, add:
```
诊断时用 SearchLKM verb="reasoning" 追溯推理链，找到 claim 的弱点和前提假设。
```

### Files NOT modified

- `main.py` (zero business logic principle)
- `src/` engine layer
- Existing tools (generate_molecules, dock_molecules, etc.)

## Scope check

- Single integration task, no sub-decomposition needed
- No architectural changes to molcraft-agent core
- Reversible: remove tool registration, .env entry, and skill directory to fully undo
