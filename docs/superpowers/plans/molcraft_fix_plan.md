# MolCraft Agent Bug Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Fix 4 critical bugs identified by Codex review that prevent proper agent execution and result submission.

**Architecture:** Incremental fixes to main.py and tools/pipeline.py, each task isolated and testable.

---

## Task 1: Fix Iteration Counting (Per-Run)

**Files:**
- Modify: `main.py:210-238` (monitor loop and count_iterations function)

**Spec:**
- At program start, record `initial_iteration_count = count_iterations()`
- In monitor, compare `count_iterations() - initial_iteration_count` against `max_iterations`
- Only trigger "已达到最大迭代次数" when NEW iterations have been completed in THIS run
- Add debug log showing iteration count delta

**Test:**
1. With existing entry in docs/iteration_log.jsonl, run `python main.py --iterations 1`
2. Verify monitor does NOT immediately stop the agent
3. Agent should proceed through all phases normally

---

## Task 2: Fix Success Validation & Stale Artifact Detection

**Files:**
- Modify: `main.py:315-360` (final metrics and packaging)

**Spec:**
- Record `run_start_time = datetime.now().timestamp()` at program start
- Before reading `result.csv`, check:
  - File exists AND
  - File mtime > run_start_time (was created/modified during THIS run)
- If no fresh result.csv:
  - Set `run_status = "failure"`
  - Log error: "No fresh pipeline output detected - run did not complete successfully"
  - Skip commit() and skip zip packaging
- Only call `StructuredLogger.commit()` and `zip_results()` if fresh result.csv exists
- result.csv path: `Path("output/result.csv")`

**Test:**
1. Delete output/result.csv, run main.py and force early termination
2. Verify run is marked as failure and no result.zip is created
3. Run successfully and verify result.zip IS created

---

## Task 3: Unify Logging - Prevent Pipeline Output Overwrite

**Files:**
- Modify: `main.py:103-140` (StructuredLogger), `tools/pipeline.py:40-80` (write_result_log)

**Spec:**
- Change `StructuredLogger` behavior: it owns the single `result.log` file
- Remove direct file writes from pipeline.py's `write_result_log()`
- Instead, pipeline logs events to stdout, which StructuredLogger captures
- OR: Pass temp log path to pipeline via env var `MOLCRAFT_LOG_PATH`
- Simpler approach: Pipeline uses same temp log file path from env
- At commit, only one file exists (no race/overwrite)

**Implementation choice:**
- In main.py, set env var `MOLCRAFT_LOG_PATH = str(tmp_log_path)` before prompt loop
- In pipeline.py, read this env var; if set, write to that path instead of default
- If env var not set (direct pipeline run), use normal `output/result.log`
- Remove the temp file overwrite pattern

**Test:**
1. Run pipeline via main.py
2. Verify both agent stdout and pipeline structured events appear in result.log
3. No data is lost or overwritten

---

## Task 4: Enforce JSONL Format Validation

**Files:**
- Modify: `main.py:135-140` (StructuredLogger), `main.py:350-360` (pre-zip validation)

**Spec:**
- In StructuredLogger.write(), wrap stdout chunks in valid JSON:
  `{"type": "stdout", "timestamp": "...", "content": "..."}`
- Write ONE complete JSON object per line (no partial lines)
- Before zip_results(), add validation step:
  - Read every line of result.log
  - Attempt json.loads() on each non-blank line
  - If any line fails parse: mark run as failure, log validation error
  - Only package valid JSONL files

**Test:**
1. Run main.py to completion
2. `python -c "import json; [json.loads(l) for l in open('output/result.log') if l.strip()]"` should succeed with no errors
