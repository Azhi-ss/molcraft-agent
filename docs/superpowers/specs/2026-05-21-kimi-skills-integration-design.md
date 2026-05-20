# Kimi Skills Integration for MolCraft Agent

**Date:** 2026-05-21
**Status:** approved

## Goal

Integrate Kimi's skill system into the molcraft-agent so the agent runtime loads domain-specific knowledge (cheminformatics, RDKit, docking strategies, retrosynthesis rules) as supplementary context without modifying `program.md`.

## Background

- `main.py` currently calls `kimi_agent_sdk.prompt()` without `skills_dir`, so no custom skills are loaded.
- kimi-cli natively supports skill discovery from `.kimi/skills/`, `.claude/skills/`, `.agents/skills/` directories.
- `skills-lock.json` already registers 4 external chemistry skills but no local SKILL.md copies exist.
- A skill is a directory containing `SKILL.md` with YAML frontmatter (`name`, `description`, optional `type`) and a markdown body. Skills are listed in the system prompt; the LLM decides when to read one.

## Directory Structure

```
molcraft-agent/
├── .kimi/
│   └── skills/
│       ├── cheminformatics/                    # external (from lockfile)
│       │   └── SKILL.md
│       ├── drug-discovery-informatics/         # external (from lockfile)
│       │   └── SKILL.md
│       ├── rdkit/                              # external (from lockfile)
│       │   └── SKILL.md
│       ├── tooluniverse-organic-chemistry/     # external (from lockfile)
│       │   └── SKILL.md
│       └── molcraft-vina-strategies/           # custom skill
│           └── SKILL.md
├── main.py                                     # modified: pass skills_dir
└── skills-lock.json                            # unchanged
```

## Changes

### 1. `main.py` — Pass `skills_dir` to `prompt()`

Add one parameter to the existing `prompt()` call:

```python
from kaos.path import KaosPath

PROJECT_ROOT = Path(__file__).parent

async for msg in prompt(
    program,
    config=llm_config,
    agent_file=AGENT_YAML,
    skills_dir=KaosPath(PROJECT_ROOT / ".kimi" / "skills"),
    yolo=True,
    thinking=args.thinking,
    max_steps_per_turn=args.max_steps,
    model=llm_model_key,
):
```

kimi-cli scans the directory, parses each `SKILL.md` frontmatter, and injects skill names + descriptions into the system prompt. The agent sees available skills and decides when to `ReadFile` a skill's full content.

### 2. External Skills — Pull from lockfile sources

Download the 4 skills registered in `skills-lock.json` from their GitHub sources:

| Skill | GitHub Source | Path |
|-------|--------------|------|
| cheminformatics | itallstartedwithaidea/agent-skills | skills/scientific-research/cheminformatics/SKILL.md |
| drug-discovery-informatics | omer-metin/skills-for-antigravity | skills/drug-discovery-informatics/SKILL.md |
| rdkit | tondevrel/scientific-agent-skills | skills/rdkit/SKILL.md |
| tooluniverse-organic-chemistry | mims-harvard/tooluniverse | skills/tooluniverse-organic-chemistry/SKILL.md |

Use `gh api` to fetch each file from GitHub raw, verify against `computedHash`, and place into `.kimi/skills/<name>/SKILL.md`.

### 3. Custom Skills — Start with one

Create `molcraft-vina-strategies` as the first custom skill, covering:
- AutoDock Vina scoring function known biases (hydrophobic overestimation, halogen underestimation)
- Docking box center/size tuning heuristics
- Common docking failure modes and countermeasures
- Target-specific knowledge (TYK2/JAK kinase hinge region characteristics)

### 4. `.gitignore` — Ensure skills are tracked

The `.kimi/skills/` directory must be tracked by Git (no gitignore entry). This ensures every experiment's skill state is reproducible and auditable.

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Directory name | `.kimi/skills/` | Natively discovered by kimi-cli, no extra config |
| Git strategy | All skills in Git | Each experiment's skill state is auditable per CLAUDE.md |
| program.md | Unchanged | Skills are supplementary knowledge, not workflow changes |
| Loading method | Explicit `skills_dir` param | More reliable than auto-discovery, self-documenting |
| Custom skill scope | Knowledge only | Does not replace tools or program.md instructions |

## Non-Goals

- Not splitting `program.md` into skills (that would change the agent workflow)
- Not building a skill manager/installer CLI (use `gh api` for one-time pull)
- Not creating flow-type skills with mermaid diagrams (standard markdown skills only)

## Verification

1. Run `main.py` and check the log for "Discovered N skill(s)" from kimi-cli
2. Confirm the system prompt includes skill names and descriptions
3. Agent can `ReadFile` a skill on demand during diagnosis phase
