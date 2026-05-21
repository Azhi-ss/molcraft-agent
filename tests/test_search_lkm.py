"""Tests for SearchLKM tool."""
import json
import os
import urllib.error
import pytest
from unittest.mock import patch


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


@pytest.mark.asyncio
async def test_missing_access_key():
    tool = _make_tool()
    params = _make_params()
    with patch.dict(os.environ, {}, clear=True):
        result = await tool(params)
    assert hasattr(result, "message")
    assert "LKM_ACCESS_KEY" in result.message or "LKM_ACCESS_KEY" in str(result.output)


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


@pytest.mark.asyncio
async def test_invalid_verb():
    tool = _make_tool()
    params = _make_params(verb="nonexistent")
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        result = await tool(params)
    assert hasattr(result, "message")


@pytest.mark.asyncio
async def test_search_missing_query():
    tool = _make_tool()
    params = _make_params(verb="search", query=None)
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        result = await tool(params)
    assert hasattr(result, "message")


@pytest.mark.asyncio
async def test_api_error_code():
    tool = _make_tool()
    params = _make_params(verb="search", query="test")
    mock_response = {
        "code": 290002,
        "message": "Invalid request parameters",
    }
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        with patch("molcraft_agent.tools._lkm_fetch_json", return_value=mock_response):
            result = await tool(params)
    output = json.loads(result.output)
    assert output["status"] == "error"
    assert output["api_code"] == 290002


@pytest.mark.asyncio
async def test_http_error():
    tool = _make_tool()
    params = _make_params(verb="search", query="test")
    err = urllib.error.HTTPError(
        url="https://example.com", code=401,
        msg="Unauthorized", hdrs=None, fp=None,
    )
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        with patch("molcraft_agent.tools._lkm_fetch_json", side_effect=err):
            result = await tool(params)
    output = json.loads(result.output)
    assert output["status"] == "error"
    assert "401" in output["error"]


@pytest.mark.asyncio
async def test_general_exception():
    tool = _make_tool()
    params = _make_params(verb="search", query="test")
    with patch.dict(os.environ, {"LKM_ACCESS_KEY": "test-key"}):
        with patch("molcraft_agent.tools._lkm_fetch_json", side_effect=ConnectionError("network down")):
            result = await tool(params)
    output = json.loads(result.output)
    assert output["status"] == "error"
    assert "network down" in output["error"]
