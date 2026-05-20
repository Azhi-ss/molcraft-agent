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
  key anchor point. The box center should be ~2-3 A from the hinge backbone toward the ATP pocket.
- **Size**: For kinase ATP sites, 20x20x20 A typically suffices. If the allosteric pocket adjacent
  to the ATP site is targeted, extend to 25x25x25 A.
- **Validation**: After setting `DOCKING_CENTER` in `src/config.py`, run a known co-crystallized
  ligand through docking. RMSD < 2.0 A from the crystal pose indicates correct box placement.

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
