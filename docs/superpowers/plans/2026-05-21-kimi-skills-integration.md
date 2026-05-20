# Kimi Skills Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire Kimi's skill system into molcraft-agent so the agent runtime loads domain knowledge skills as supplementary context.

**Architecture:** Create `.kimi/skills/` with 4 external skills (pulled from GitHub via `skills-lock.json` sources) + 1 custom skill (`molcraft-vina-strategies`). Modify `main.py` to pass `skills_dir=KaosPath(...)` to the existing `prompt()` call. kimi-cli auto-discovers SKILL.md files and injects them into the system prompt.

**Tech Stack:** Python 3.12, Kimi Agent SDK, kimi-cli, `gh` CLI, KaosPath

---

### Task 1: Create `.kimi/skills/` directory structure

**Files:**
- Create: `.kimi/skills/cheminformatics/SKILL.md`
- Create: `.kimi/skills/drug-discovery-informatics/SKILL.md`
- Create: `.kimi/skills/rdkit/SKILL.md`
- Create: `.kimi/skills/tooluniverse-organic-chemistry/SKILL.md`

- [ ] **Step 1: Create directories**

```bash
mkdir -p .kimi/skills/cheminformatics
mkdir -p .kimi/skills/drug-discovery-informatics
mkdir -p .kimi/skills/rdkit
mkdir -p .kimi/skills/tooluniverse-organic-chemistry
```

- [ ] **Step 2: Fetch cheminformatics skill from GitHub**

```bash
gh api repos/itallstartedwithaidea/agent-skills/contents/skills/scientific-research/cheminformatics/SKILL.md \
  --jq '.content' | base64 -d > .kimi/skills/cheminformatics/SKILL.md
```

- [ ] **Step 3: Fetch drug-discovery-informatics skill from GitHub**

```bash
gh api repos/omer-metin/skills-for-antigravity/contents/skills/drug-discovery-informatics/SKILL.md \
  --jq '.content' | base64 -d > .kimi/skills/drug-discovery-informatics/SKILL.md
```

- [ ] **Step 4: Fetch rdkit skill from GitHub**

```bash
gh api repos/tondevrel/scientific-agent-skills/contents/skills/rdkit/SKILL.md \
  --jq '.content' | base64 -d > .kimi/skills/rdkit/SKILL.md
```

- [ ] **Step 5: Fetch tooluniverse-organic-chemistry skill from GitHub**

```bash
gh api repos/mims-harvard/tooluniverse/contents/skills/tooluniverse-organic-chemistry/SKILL.md \
  --jq '.content' | base64 -d > .kimi/skills/tooluniverse-organic-chemistry/SKILL.md
```

- [ ] **Step 6: Verify hashes match skills-lock.json**

```bash
echo "cheminformatics: $(sha256sum .kimi/skills/cheminformatics/SKILL.md | cut -d' ' -f1)"
echo "expected: 9f96c04e34871585f1dbbb3ad4e41d32b709a4318aa350fcc9ba52d70c4fbc7e"
echo ""
echo "drug-discovery-informatics: $(sha256sum .kimi/skills/drug-discovery-informatics/SKILL.md | cut -d' ' -f1)"
echo "expected: bdca6426b0833de85511d17281cd385b097908fb8567b30daff799ff2de2ae0b"
echo ""
echo "rdkit: $(sha256sum .kimi/skills/rdkit/SKILL.md | cut -d' ' -f1)"
echo "expected: 160e9a8be4feaee922de321819d0eccd6df165cec19501bc6cdb0746fa77ed1e"
echo ""
echo "tooluniverse-organic-chemistry: $(sha256sum .kimi/skills/tooluniverse-organic-chemistry/SKILL.md | cut -d' ' -f1)"
echo "expected: 32c2b681dbae220ce7256519843c0c72c71edbabddbf37d9a7b1072d2052b0ec"
```

- [ ] **Step 7: Commit external skills**

```bash
git add .kimi/skills/cheminformatics/SKILL.md \
        .kimi/skills/drug-discovery-informatics/SKILL.md \
        .kimi/skills/rdkit/SKILL.md \
        .kimi/skills/tooluniverse-organic-chemistry/SKILL.md
git commit -m "feat: add 4 external chemistry skills to .kimi/skills/

cheminformatics, drug-discovery-informatics, rdkit, tooluniverse-organic-chemistry
pulled from GitHub sources per skills-lock.json hashes.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 2: Wire `skills_dir` into `main.py`

**Files:**
- Modify: `main.py:17,301-308`

- [ ] **Step 1: Add KaosPath import to main.py**

At line 27 (after `from pathlib import Path`), add:

```python
from kaos.path import KaosPath
```

- [ ] **Step 2: Add PROJECT_ROOT constant**

At line 32 (after imports, before `AGENT_YAML`), add:

```python
PROJECT_ROOT = Path(__file__).parent
```

- [ ] **Step 3: Add skills_dir to the prompt() call**

At line 304 (inside the `prompt()` kwargs), insert after `agent_file=AGENT_YAML,`:

```python
            skills_dir=KaosPath(PROJECT_ROOT / ".kimi" / "skills"),
```

The resulting `prompt()` call should read:

```python
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

- [ ] **Step 4: Verify syntax — dry-run import**

```bash
.venv/bin/python -c "import ast; ast.parse(open('main.py').read()); print('Syntax OK')"
```
Expected: `Syntax OK`

- [ ] **Step 5: Commit**

```bash
git add main.py
git commit -m "feat: pass skills_dir to prompt() for Kimi skill loading

Adds KaosPath import and PROJECT_ROOT constant, wires .kimi/skills/
as the skills directory so the agent runtime discovers and lists
available skills in the system prompt.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 3: Create custom skill `molcraft-vina-strategies`

**Files:**
- Create: `.kimi/skills/molcraft-vina-strategies/SKILL.md`

- [ ] **Step 1: Create skill directory**

```bash
mkdir -p .kimi/skills/molcraft-vina-strategies
```

- [ ] **Step 2: Write SKILL.md**

```markdown
---
name: molcraft-vina-strategies
description: AutoDock Vina docking optimization knowledge for TYK2/JAK kinase targets. Covers scoring function biases, docking box tuning, common failure modes, and kinase-specific binding pocket characteristics.
---

# Vina Docking Strategies for Kinase Targets

## Scoring Function Known Biases

- **Hydrophobic overestimation**: Vina over-rewards burial of nonpolar surface area. Large aromatic
  substituents in solvent-exposed regions may score well but bind poorly in reality.
- **Halogen underestimation**: Halogen bonds (C-X...O=C) contribute ~0.5-2 kcal/mol but Vina
  lacks explicit sigma-hole terms. Halogenated fragments near backbone carbonyls are under-scored.
- **Desolvation penalty weak**: Vina's implicit solvent model underestimates the cost of
  desolvating polar groups. Highly polar molecules may appear better than they are.
- **Entropic penalty flat**: The ligand entropy penalty is approximated as proportional to the
  number of rotatable bonds (~0.3 kcal/mol per bond). This is a rough estimate.

## Docking Box Tuning Heuristics

- **Center**: Must be at the geometric center of the binding site, NOT the protein center of mass.
  For TYK2 JH1 domain, the hinge region (residues 950-960, backbone NH of Met953, Glu955) is the
  key anchor point. The box center should be ~2-3 Å from the hinge backbone toward the ATP pocket.
- **Size**: For kinase ATP sites, 20×20×20 Å typically suffices. If the allosteric pocket adjacent
  to the ATP site is targeted, extend to 25×25×25 Å.
- **Validation**: After setting `DOCKING_CENTER` in `src/config.py`, run a known co-crystallized
  ligand through docking. RMSD < 2.0 Å from the crystal pose indicates correct box placement.

## Common Docking Failure Modes

| Symptom | Likely Cause | Fix |
|---------|-------------|-----|
| All scores > -5 kcal/mol | Box center off-target | Verify against crystal structure pocket residues |
| Zero successful docks | Molecule too large for box | Increase box size or filter by MW < 500 |
| All scores clustered at -4 to -6 | Exhaustiveness too low | Increase `--exhaustiveness` from 8 to 16 |
| Good scores but bad poses | Scoring function deceived by hydrophobics | Cross-validate with pose quality check |

## Kinase Hinge Region Characteristics

- **Hinge motif**: The backbone NH of the hinge residue (Met953 in TYK2) is a conserved hydrogen
  bond donor. A heterocyclic nitrogen (pyridine, pyrimidine, pyrazole) that accepts this H-bond
  is the hallmark of type I kinase inhibitors.
- **Gatekeeper residue**: TYK2 has Thr972 as gatekeeper. Bulkier gatekeepers (Phe, Met) restrict
  access to the back pocket; smaller ones (Thr) allow type II inhibitors to extend past the gate.
- **DFG motif**: Asp-Phe-Gly at the activation loop N-terminus. DFG-in = ATP-competitive (type I),
  DFG-out = allosteric (type II). Vina cannot distinguish these states -- always use the DFG-in
  conformation for docking unless explicitly targeting the allosteric pocket.
- **Ribose pocket**: The conserved Asp1041 interacts with the ribose hydroxyls of ATP. Polar
  fragments (alcohol, amide) near this residue can improve selectivity.

## Recommended Docking Protocol for This Project

1. Verify `DOCKING_CENTER` in `src/config.py` with `identify_target()`
2. Run a small batch (n=20) as smoke test: check docking success rate > 80%
3. If success rate < 80%, check box size and molecule MW distribution
4. After each generator change, re-run smoke test to confirm docking still works
5. Cross-reference top 10% binding energies with QED: if best BE comes from QED < 0.3 molecules,
   the scoring function may be overweighting unrealistic features
```

- [ ] **Step 3: Validate SKILL.md frontmatter is parseable**

```bash
.venv/bin/python -c "
from kimi_cli.skill import parse_skill_text
from kaos.path import KaosPath
skill = parse_skill_text(open('.kimi/skills/molcraft-vina-strategies/SKILL.md').read(), dir_path=KaosPath('.kimi/skills/molcraft-vina-strategies'))
print(f'name={skill.name}')
print(f'description={skill.description[:50]}...')
print(f'type={skill.type}')
"
```
Expected: prints `name=molcraft-vina-strategies`, description snippet, `type=standard`

- [ ] **Step 4: Commit**

```bash
git add .kimi/skills/molcraft-vina-strategies/SKILL.md
git commit -m "feat: add molcraft-vina-strategies custom skill

Covers Vina scoring biases, docking box tuning, kinase hinge
characteristics, and common failure modes for TYK2/JAK targets.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 4: End-to-end verification

- [ ] **Step 1: Confirm skill discovery at agent startup**

Run a minimal agent invocation (1 step, dry-run mode if available) and check the log for skill discovery:

```bash
.venv/bin/python -c "
import asyncio
from pathlib import Path
from kaos.path import KaosPath
from kimi_agent_sdk import prompt

async def test():
    skills_dir = KaosPath(Path('.kimi/skills').resolve())
    count = 0
    async for msg in prompt(
        'echo skills loaded',
        skills_dir=skills_dir,
        yolo=True,
        max_steps_per_turn=1,
        final_message_only=True,
    ):
        count += 1
    print(f'Messages received: {count}')

asyncio.run(test())
" 2>&1 | head -5
```

- [ ] **Step 2: Verify .gitignore does not exclude .kimi/**

```bash
git check-ignore -v .kimi/skills/cheminformatics/SKILL.md
```
Expected: no output (file is NOT ignored)

- [ ] **Step 3: Verify git status clean after all commits**

```bash
git status
```
Expected: `nothing to commit, working tree clean`
