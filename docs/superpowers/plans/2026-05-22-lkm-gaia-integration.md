# LKM + Gaia Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Integrate Bohrium LKM scientific knowledge graph search and Gaia formalize into molcraft-agent, giving the Agent structured literature search and knowledge-package creation capabilities.

**Architecture:** New `SearchLKM` Python tool for high-frequency LKM API calls (search/reasoning/papers-graph). Gaia formalize skills copied as Kimi skill directories under `.kimi/skills/lkm/`. Config changes to `.env`, `agent.yaml`, and `program.md`.

**Tech Stack:** Python 3, urllib.request (stdlib), Pydantic, Kimi Agent SDK (CallableTool2), Gaia CLI

---

## File Structure

| File | Action | Responsibility |
|------|--------|---------------|
| `molcraft_agent/tools.py` | Modify | Add `SearchLKM` tool class + params model |
| `.env` | Modify | Add `LKM_ACCESS_KEY` |
| `agent.yaml` | Modify | Register `SearchLKM` tool |
| `program.md` | Modify | Add LKM search priority in stages 1 & 2 |
| `.kimi/skills/lkm/orchestrator/SKILL.md` | Create (copy) | LKM/Gaia skill router |
| `.kimi/skills/lkm/orchestrator/references/` | Create (copy) | SOP references |
| `.kimi/skills/lkm/lkm-search/SKILL.md` | Create (copy) | LKM API contract skill |
| `.kimi/skills/lkm/lkm-search/scripts/lkm_search.py` | Create (copy) | CLI helper for low-freq endpoints |
| `.kimi/skills/lkm/lkm-search/references/api-contract.md` | Create (copy) | API schema reference |
| `.kimi/skills/lkm/lkm-search/.gitignore` | Create (copy) | Skill version marker ignore |
| `.kimi/skills/lkm/lkm-explorer/SKILL.md` | Create (copy) | LKM→Gaia workflow skill |
| `.kimi/skills/lkm/lkm-explorer/references/` | Create (copy) | Explorer step references |
| `.kimi/skills/lkm/formalize/SKILL.md` | Create (copy) | Paper→Gaia workflow skill |
| `.kimi/skills/lkm/formalize/references/` | Create (copy) | Formalize phase references |
| `tests/test_search_lkm.py` | Create | Unit tests for SearchLKM tool |

---

### Task 1: Add SearchLKM Tool

**Files:**
- Modify: `molcraft_agent/tools.py` (append after `DiffusionGenerate` class, ~line 1180)
- Create: `tests/test_search_lkm.py`

- [ ] **Step 1: Write the test file**

Create `tests/test_search_lkm.py`:

```python
"""Tests for SearchLKM tool."""
import json
import os
import pytest
from unittest.mock import patch, MagicMock


# -- Helper to instantiate the tool --

def _make_tool():
    from molcraft_agent.tools import SearchLKM
    return SearchLKM()


def _make_params(**overrides):
    from molcraft_agent.tools import SearchLKMParams
    defaults = {
        "verb": "search",
        "query": "kinase inhibitor",
        "claim_id": None,
        "doi": None,
        "title": None,
        "top_k": 20,
        "reasoning_only": False,
    }
    defaults.update(overrides)
    return SearchLKMParams(**defaults)


# -- Missing API key --

@pytest.mark.asyncio
async def test_missing_access_key():
    tool = _make_tool()
    params = _make_params()
    with patch.dict(os.environ, {}, clear=True):
        result = await tool(params)
    # Should return ToolError
    assert hasattr(result, "message")
    assert "LKM_ACCESS_KEY" in result.message or "LKM_ACCESS_KEY" in str(result.output)


# -- Search verb --

@pytest.mark.asyncio
async def test_search_verb_success():
    tool = _make_tool()
    params = _make_params(verb="search", query="perovskite stability")
    mock_response = {
        "code": 0,
        "data": {
            "variables": [{"id": "gcn_001", "name": "test claim"}],
            "total": 1,
        },
    }
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        with patch("molcraft_agent.tools._lkm_fetch_json", return_value=mock_response):
            result = await tool(params)
    output = json.loads(result.output)
    assert output["status"] == "success"
    assert output["verb"] == "search"
    assert "data" in output


# -- Reasoning verb --

@pytest.mark.asyncio
async def test_reasoning_verb_success():
    tool = _make_tool()
    params = _make_params(verb="reasoning", claim_id="gcn_abc123")
    mock_response = {
        "code": 0,
        "data": {
            "chains": [{"id": "chain_1"}],
            "total_chains": 1,
        },
    }
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        with patch("molcraft_agent.tools._lkm_fetch_json", return_value=mock_response):
            result = await tool(params)
    output = json.loads(result.output)
    assert output["status"] == "success"
    assert output["verb"] == "reasoning"


# -- Papers graph verb --

@pytest.mark.asyncio
async def test_papers_graph_verb_success():
    tool = _make_tool()
    params = _make_params(verb="papers_graph", doi="10.1038/s41586-023-06408-7")
    mock_response = {
        "code": 0,
        "data": {
            "paper": {"title": "Test Paper"},
            "variables": [],
        },
    }
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        with patch("molcraft_agent.tools._lkm_fetch_json", return_value=mock_response):
            result = await tool(params)
    output = json.loads(result.output)
    assert output["status"] == "success"
    assert output["verb"] == "papers_graph"


# -- Invalid verb --

@pytest.mark.asyncio
async def test_invalid_verb():
    tool = _make_tool()
    from molcraft_agent.tools import SearchLKMParams
    # Pydantic should accept any string; validation is in __call__
    params = _make_params(verb="nonexistent")
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        result = await tool(params)
    assert hasattr(result, "message")  # ToolError


# -- Missing required param for verb --

@pytest.mark.asyncio
async def test_search_missing_query():
    tool = _make_tool()
    params = _make_params(verb="search", query=None)
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        result = await tool(params)
    assert hasattr(result, "message")  # ToolError
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /home/dministrator/Lab/clones/molcraft-agent && python -m pytest tests/test_search_lkm.py -v 2>&1 | tail -20`
Expected: ImportError or ModuleNotFoundError — `SearchLKM` and `SearchLKMParams` not defined yet.

- [ ] **Step 3: Add SearchLKMParams model and SearchLKM tool to tools.py**

Append the following to the end of `molcraft_agent/tools.py` (after line 1180):

```python
# ── LKM Knowledge Graph Search ──

LKM_BASE_URL = "https://open.bohrium.com/openapi/v1/lkm"


def _lkm_fetch_json(
    url: str,
    *,
    method: str = "GET",
    headers: dict[str, str] | None = None,
    body: Any | None = None,
) -> Any:
    """HTTP request to Bohrium LKM API. Returns parsed JSON."""
    data: bytes | None = None
    if body is not None:
        data = json.dumps(body).encode("utf-8")
    req = urllib.request.Request(
        url, data=data, headers=headers or {}, method=method,
    )
    with urllib.request.urlopen(req) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    return json.loads(text)


class SearchLKMParams(BaseModel):
    verb: str = Field(
        description="子命令: search(搜索科学claims), reasoning(追溯推理链), papers_graph(论文知识图谱)",
    )
    query: str | None = Field(
        default=None,
        description="搜索查询（verb=search 时必填）",
    )
    claim_id: str | None = Field(
        default=None,
        description="Claim GCN ID（verb=reasoning 时必填）",
    )
    doi: str | None = Field(
        default=None,
        description="论文 DOI（verb=papers_graph 时使用）",
    )
    title: str | None = Field(
        default=None,
        description="论文标题（verb=papers_graph 时使用）",
    )
    top_k: int = Field(
        default=20,
        description="搜索返回最大结果数",
    )
    reasoning_only: bool = Field(
        default=False,
        description="仅返回有推理链的 claims",
    )


class SearchLKM(CallableTool2):
    name: str = "search_lkm"
    description: str = (
        "搜索 Bohrium LKM 科学知识图谱。三个子命令："
        "search — 按主题搜索科学 claims 和研究问题；"
        "reasoning — 追溯某条 claim 的推理链（前提、推导步骤、弱点）；"
        "papers_graph — 获取论文的完整知识图谱。"
        "文献解析阶段优先使用此工具，SearchWeb 作为兜底。"
    )
    params: type[BaseModel] = SearchLKMParams

    async def __call__(self, params: SearchLKMParams) -> ToolReturnValue:
        access_key = os.environ.get("LKM_ACCESS_KEY")
        if not access_key:
            return ToolError(
                output="",
                message="LKM_ACCESS_KEY 未设置。请在 .env 中添加 LKM_ACCESS_KEY=<key>，或 export LKM_ACCESS_KEY=<key>",
                brief="LKM 未授权",
            )

        try:
            if params.verb == "search":
                result_data = await asyncio.to_thread(
                    self._search, access_key, params,
                )
            elif params.verb == "reasoning":
                result_data = await asyncio.to_thread(
                    self._reasoning, access_key, params,
                )
            elif params.verb == "papers_graph":
                result_data = await asyncio.to_thread(
                    self._papers_graph, access_key, params,
                )
            else:
                return ToolError(
                    output="",
                    message=f"未知的 verb: {params.verb}。支持: search, reasoning, papers_graph",
                    brief="无效子命令",
                )

            output = {
                "status": "success",
                "verb": params.verb,
                "data": result_data,
            }
            append_experiment(
                tool="search_lkm",
                round_num=get_latest_round(),
                params={"verb": params.verb, "query": params.query, "claim_id": params.claim_id},
                result={"code": result_data.get("code"), "total": result_data.get("data", {}).get("total") if isinstance(result_data.get("data"), dict) else None},
            )
            return ToolOk(output=json.dumps(output, ensure_ascii=False))

        except urllib.error.HTTPError as exc:
            body_text = exc.read().decode("utf-8", errors="replace")[:500]
            return ToolError(
                output=json.dumps({"status": "error", "verb": params.verb, "http_status": exc.code, "detail": body_text}),
                message=f"LKM API HTTP {exc.code}",
                brief="LKM 请求失败",
            )
        except Exception as exc:
            return ToolError(
                output=json.dumps({"status": "error", "verb": params.verb, "error": str(exc)}),
                message=str(exc),
                brief="LKM 搜索失败",
            )

    @staticmethod
    def _search(access_key: str, params: SearchLKMParams) -> dict:
        if not params.query:
            raise ValueError("verb=search 需要 query 参数")
        body: dict[str, Any] = {
            "query": params.query,
            "top_k": params.top_k,
        }
        if params.reasoning_only:
            body["reasoning_only"] = True
        return _lkm_fetch_json(
            f"{LKM_BASE_URL}/search",
            method="POST",
            headers={
                "accessKey": access_key,
                "content-type": "application/json",
            },
            body=body,
        )

    @staticmethod
    def _reasoning(access_key: str, params: SearchLKMParams) -> dict:
        if not params.claim_id:
            raise ValueError("verb=reasoning 需要 claim_id 参数")
        claim_id = urllib.parse.quote(params.claim_id, safe="")
        url = f"{LKM_BASE_URL}/claims/{claim_id}/reasoning?max_chains=10&sort_by=comprehensive"
        return _lkm_fetch_json(
            url,
            headers={"accessKey": access_key},
        )

    @staticmethod
    def _papers_graph(access_key: str, params: SearchLKMParams) -> dict:
        body: dict[str, Any] = {}
        if params.doi:
            body["doi"] = params.doi
        if params.title:
            body["title"] = params.title
        if not body:
            raise ValueError("verb=papers_graph 需要 doi 或 title 参数")
        return _lkm_fetch_json(
            f"{LKM_BASE_URL}/papers/graph",
            method="POST",
            headers={
                "accessKey": access_key,
                "content-type": "application/json",
            },
            body=body,
        )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/dministrator/Lab/clones/molcraft-agent && python -m pytest tests/test_search_lkm.py -v 2>&1 | tail -30`
Expected: All 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add molcraft_agent/tools.py tests/test_search_lkm.py
git commit -m "feat: add SearchLKM tool for Bohrium LKM knowledge graph"
```

---

### Task 2: Register SearchLKM in agent.yaml

**Files:**
- Modify: `agent.yaml` (after line 29, the DiffusionGenerate entry)

- [ ] **Step 1: Add tool registration**

In `agent.yaml`, after the line `- "molcraft_agent.tools:DiffusionGenerate"` (line 29), add:

```yaml
    # LKM 科学知识图谱搜索
    - "molcraft_agent.tools:SearchLKM"
```

- [ ] **Step 2: Verify YAML is valid**

Run: `cd /home/dministrator/Lab/clones/molcraft-agent && python -c "import yaml; yaml.safe_load(open('agent.yaml'))"`
Expected: No error.

- [ ] **Step 3: Commit**

```bash
git add agent.yaml
git commit -m "feat: register SearchLKM tool in agent.yaml"
```

---

### Task 3: Add LKM_ACCESS_KEY to .env

**Files:**
- Modify: `.env`

- [ ] **Step 1: Append key to .env**

Append to the end of `.env`:

```
LKM_ACCESS_KEY="4fc2fbe980fe452bbb0c8d09d80f056a"
```

- [ ] **Step 2: Verify .env loads**

Run: `cd /home/dministrator/Lab/clones/molcraft-agent && python -c "from dotenv import load_dotenv; load_dotenv(); import os; print('LKM_ACCESS_KEY' in os.environ)"`
Expected: `True`

- [ ] **Step 3: Commit**

```bash
git add .env
git commit -m "feat: add LKM_ACCESS_KEY to .env"
```

---

### Task 4: Copy Gaia-LKM Skills

**Files:**
- Create: `.kimi/skills/lkm/` (multiple subdirectories and files)

- [ ] **Step 1: Copy orchestrator skill**

```bash
mkdir -p /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/orchestrator/references
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/orchestrator/SKILL.md /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/orchestrator/
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/orchestrator/references/lkm-explorer-sop.md /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/orchestrator/references/
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/orchestrator/references/audited-delegation.md /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/orchestrator/references/
```

- [ ] **Step 2: Copy lkm-search skill**

```bash
mkdir -p /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/lkm-search/scripts
mkdir -p /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/lkm-search/references
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/lkm-search/SKILL.md /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/lkm-search/
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/lkm-search/.gitignore /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/lkm-search/
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/lkm-search/scripts/lkm_search.py /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/lkm-search/scripts/
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/lkm-search/references/api-contract.md /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/lkm-search/references/
```

- [ ] **Step 3: Copy lkm-explorer skill**

```bash
mkdir -p /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/lkm-explorer/references
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/lkm-explorer/SKILL.md /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/lkm-explorer/
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/lkm-explorer/references/*.md /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/lkm-explorer/references/
```

- [ ] **Step 4: Copy formalize skill**

```bash
mkdir -p /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/formalize/references
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/formalize/SKILL.md /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/formalize/
cp /home/dministrator/Lab/clones/gaia-lkm-skills/skills/formalize/references/*.md /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/formalize/references/
```

- [ ] **Step 5: Verify skill files are in place**

Run: `find /home/dministrator/Lab/clones/molcraft-agent/.kimi/skills/lkm/ -type f | sort`
Expected: 16 files (1 orchestrator SKILL + 2 references, 1 lkm-search SKILL + 1 gitignore + 1 script + 1 reference, 1 lkm-explorer SKILL + 7 references, 1 formalize SKILL + 4 references = 20 files)

- [ ] **Step 6: Commit**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent
git add .kimi/skills/lkm/
git commit -m "feat: add Gaia-LKM skills (orchestrator, lkm-search, lkm-explorer, formalize)"
```

---

### Task 5: Update program.md with LKM Search Priority

**Files:**
- Modify: `program.md`

- [ ] **Step 1: Add LKM priority in 阶段一 (line 52, after "首次运行" paragraph)**

Find the paragraph starting with `**首次运行**：读 papers/ 中的论文` (line 52). Insert the following block right before that paragraph:

```
**文献检索优先级**：
1. `SearchLKM` — 搜索 Bohrium LKM 科学知识图谱，获取结构化 claims 和推理链（优先使用）
2. `SearchWeb` — LKM 无结果时兜底，搜通用文献和博客
3. 当需要将 LKM 证据图谱化时，读取 `.kimi/skills/lkm/orchestrator/SKILL.md` 按 SOP 走 `lkm-explorer` 或 `formalize`

```

- [ ] **Step 2: Add reasoning guidance in 阶段二 (line 59, after "禁重检查")**

Find the line `0. **禁重检查（每轮必做）**` (line 60). After the禁重检查block, and before `1. **搜索外部新知识`, add a new bullet after 搜索外部新知识 that references LKM reasoning. Replace the existing bullet `1. **搜索外部新知识（每轮必做）**` content to include LKM:

In the line starting `1. **搜索外部新知识（每轮必做）**：搜索当前靶点/瓶颈的最新进展。`, change it to:

```
1. **搜索外部新知识（每轮必做）**：优先用 `SearchLKM verb="search"` 搜索 LKM 知识图谱中当前靶点/瓶颈的最新进展，提取至少 1 个与知识库对比过的改进方向。用 `SearchLKM verb="reasoning"` 追溯推理链，找到 claim 的弱点和前提假设。LKM 无结果时再用 SearchWeb 兜底。
```

- [ ] **Step 3: Verify program.md is well-formed**

Run: `cd /home/dministrator/Lab/clones/molcraft-agent && python -c "open('program.md').read(); print('OK')"`
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add program.md
git commit -m "feat: add LKM search priority to program.md stages 1 & 2"
```

---

### Task 6: Smoke Test — Verify SearchLKM Tool Works Against Live API

**Files:**
- No file changes (validation only)

- [ ] **Step 1: Run a real search query**

```bash
cd /home/dministrator/Lab/clones/molcraft-agent
source .venv/bin/activate
python -c "
from dotenv import load_dotenv; load_dotenv()
import os, json
from molcraft_agent.tools import SearchLKM, SearchLKMParams
import asyncio

async def test():
    tool = SearchLKM()
    params = SearchLKMParams(verb='search', query='kinase inhibitor binding', top_k=5)
    result = await tool(params)
    print(type(result).__name__)
    if hasattr(result, 'output'):
        data = json.loads(result.output)
        print(f'status: {data.get(\"status\")}')
        if data.get('data', {}).get('data', {}).get('variables'):
            print(f'found {len(data[\"data\"][\"data\"][\"variables\"])} results')
        else:
            print(json.dumps(data, indent=2, ensure_ascii=False)[:500])

asyncio.run(test())
"
```

Expected: `ToolOk`, `status: success`, and some search results from LKM.

- [ ] **Step 2: If API returns error, diagnose and fix**

Common issues:
- HTTP 401: accessKey wrong → check `.env` LKM_ACCESS_KEY value
- HTTP 400: query format → adjust default params
- Timeout: network issue → retry once

---

## Self-Review Checklist

1. **Spec coverage**:
   - SearchLKM tool with 3 verbs → Task 1 ✓
   - Gaia formalize skills copied → Task 4 ✓
   - .env LKM_ACCESS_KEY → Task 3 ✓
   - agent.yaml registration → Task 2 ✓
   - program.md search priority → Task 5 ✓
   - Smoke test → Task 6 ✓

2. **Placeholder scan**: No TBD/TODO/implement-later found. All code blocks are complete.

3. **Type consistency**: `SearchLKMParams` fields match between test and implementation. `_lkm_fetch_json` signature used consistently in both tests (mock target) and implementation.
