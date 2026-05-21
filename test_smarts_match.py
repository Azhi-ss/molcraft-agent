"""
Diagnose SMARTS mis-match for molecules #2, #3, #4 in result.csv
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# The target molecules from result.csv (full molecules)
targets = {
    "#2: Cc1ccc(-c2cncc3ccccc23)cc1": "Cc1ccc(-c2cncc3ccccc23)cc1",
    "#3: FCc1ncc2cccc(-c3ccccc3)c2n1": "FCc1ncc2cccc(-c3ccccc3)c2n1",
    "#4: c1ccc(-c2cccc3ncncc23)cc1": "c1ccc(-c2cccc3ncncc23)cc1",
}

# The intermediates from the first retro step (product side of step 1)
intermediates = {
    "#2 first step product": "OB(O)c1cncc2ccccc12",
    "#3 first step product": "FCc1ncc2cccc(Br)c2n1",
    "#4 first step product": "OB(O)c1cccc2ncncc12",
}

# All SMARTS patterns from RETRO_RULES that produce C=O
retro_rules_with_co = [
    # line 239: isoquinoline → Pictet-Spengler
    ("H015: isoquinoline >> NCCc1ccccc1.C=O", 
     "c1ccc2c(c1)ccnc2", "c1ccc2c(c1)ccnc2>>NCCc1ccccc1.C=O"),
    # line 237: quinoline → Friedländer
    ("H012: quinoline >> Nc1ccccc1C=O.CC=O",
     "c1ccc2c(c1)cccn2", "c1ccc2c(c1)cccn2>>Nc1ccccc1C=O.CC=O"),
    # line 243: quinazoline → anthranilamide + formaldehyde
    ("H012: quinazoline >> Nc1ccccc1C(=O)N.C=O",
     "c1ccc2c(c1)ncnc2", "c1ccc2c(c1)ncnc2>>Nc1ccccc1C(=O)N.C=O"),
    # line 344: THIQ >> phenethylamine + C=O
    ("H015: THIQ >> NCCc1ccccc1.C=O",
     "c1ccc2c(c1)CCNC2", "c1ccc2c(c1)CCNC2>>NCCc1ccccc1.C=O"),
    # line 493: indole >> phenylhydrazine + C=O
    ("H020: indole >> NNc1ccccc1.C=O",
     "c1ccc2[nH]ccc2c1", "c1ccc2[nH]ccc2c1>>NNc1ccccc1.C=O"),
]

print("=" * 100)
print("TEST 1: Which SMARTS patterns match the intermediates (first retro step products)?")
print("=" * 100)

for name, smi in intermediates.items():
    print(f"\n--- {name} (SMILES: {smi}) ---")
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  ERROR: Cannot parse SMILES")
        continue
    print(f"  Heavy atoms: {mol.GetNumHeavyAtoms()}, Total atoms: {mol.GetNumAtoms()}")
    
    # Show the ring info
    ri = mol.GetRingInfo()
    print(f"  Rings: {ri.NumRings()}")
    
    for rule_name, smarts_pat, retro_smarts in retro_rules_with_co:
        patt = Chem.MolFromSmarts(smarts_pat)
        if patt is None:
            print(f"  Rule [{rule_name}]: SMARTS parse error")
            continue
        match = mol.HasSubstructMatch(patt)
        if match:
            matches = mol.GetSubstructMatches(patt)
            print(f"  ✓ MATCHED by [{rule_name}]")
            print(f"    SMARTS pattern: {smarts_pat}")
            
            # Try running the retro reaction
            rxn = AllChem.ReactionFromSmarts(retro_smarts)
            if rxn:
                ps = rxn.RunReactants((mol,))
                if ps and ps[0]:
                    product_smiles = [Chem.MolToSmiles(p) for p in ps[0]]
                    print(f"    Reactants produced: {product_smiles}")
                    
                    # Check atoms
                    for p_smi in product_smiles:
                        p_mol = Chem.MolFromSmiles(p_smi)
                        if p_mol:
                            print(f"      {p_smi}: {p_mol.GetNumHeavyAtoms()} heavy atoms")
                    
                    # Check if atoms sum up reasonably
                    r_heavy = sum(Chem.MolFromSmiles(p).GetNumHeavyAtoms() for p in product_smiles if Chem.MolFromSmiles(p))
                    t_heavy = mol.GetNumHeavyAtoms()
                    ratio = t_heavy / r_heavy if r_heavy > 0 else 0
                    print(f"    Target heavy: {t_heavy}, Reactant heavy: {r_heavy}, Ratio: {ratio:.2f}")
                    if ratio > 1.3:
                        print(f"    ✗ ATOM IMBALANCE: target has {ratio:.2f}x heavy atoms of reactants")

print("\n" + "=" * 100)
print("TEST 2: Do the intermediates contain the heterocyclic cores?")
print("=" * 100)

# The core patterns
core_patterns = {
    "isoquinoline (c1ccc2c(c1)ccnc2)": "c1ccc2c(c1)ccnc2",
    "quinoline (c1ccc2c(c1)cccn2)": "c1ccc2c(c1)cccn2",
    "quinazoline (c1ccc2c(c1)ncnc2)": "c1ccc2c(c1)ncnc2",
    "indole (c1ccc2[nH]ccc2c1)": "c1ccc2[nH]ccc2c1",
}

for iname, ismi in intermediates.items():
    imol = Chem.MolFromSmiles(ismi)
    print(f"\n{iname}: {ismi}")
    for cname, csmarts in core_patterns.items():
        cpatt = Chem.MolFromSmarts(csmarts)
        if cpatt and imol and imol.HasSubstructMatch(cpatt):
            print(f"  CONTAINS core: {cname}")

print("\n" + "=" * 100)
print("TEST 3: What does quinazoline SMARTS really match?")
print("=" * 100)

# The quinazoline core
quinazoline_smiles = "c1ncnc2ccccc12"
quinazoline_mol = Chem.MolFromSmiles(quinazoline_smiles)
print(f"Quinazoline canonical SMILES: {Chem.MolToSmiles(quinazoline_mol)}")

# Test what intermediates match the quinazoline SMARTS
qz_patt = Chem.MolFromSmarts("c1ccc2c(c1)ncnc2")
for iname, ismi in intermediates.items():
    imol = Chem.MolFromSmiles(ismi)
    if imol and qz_patt:
        matches = imol.GetSubstructMatches(qz_patt)
        print(f"\n{iname}: {ismi}")
        print(f"  Matches quinazoline core: {imol.HasSubstructMatch(qz_patt)}")
        if matches:
            print(f"  Match positions: {matches}")
            # Show atoms in the match
            for match in matches:
                for atom_idx in match:
                    atom = imol.GetAtomWithIdx(atom_idx)
                    print(f"    Atom {atom_idx}: {atom.GetSymbol()} (degree={atom.GetDegree()}, neighbors={[n.GetIdx() for n in atom.GetNeighbors()]})")

# Also check: does FCc1ncc2cccc(Br)c2n1 have a naphthyridine-like structure?
print("\n" + "=" * 100)
print("TEST 4: Detailed structure of problem molecules")
print("=" * 100)

for name, smi in list(intermediates.items()):
    mol = Chem.MolFromSmiles(smi)
    print(f"\n{name}: {smi}")
    print(f"  Canonical: {Chem.MolToSmiles(mol)}")
    print(f"  Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
    print(f"  Heavy atoms: {mol.GetNumHeavyAtoms()}")
    
    # List all atoms with their symbols
    for i, atom in enumerate(mol.GetAtoms()):
        print(f"    Atom {i}: {atom.GetSymbol()} (ring={atom.IsInRing()}, aromatic={atom.GetIsAromatic()})")

# Check molecule #2 intermediate: is it REALLY isoquinoline?
print("\n" + "=" * 100)
print("TEST 5: Molecule #2 intermediate - is isoquinoline match correct?")
print("=" * 100)

smi2 = "OB(O)c1cncc2ccccc12"
mol2 = Chem.MolFromSmiles(smi2)
print(f"Molecule: {smi2}")
print(f"Canonical: {Chem.MolToSmiles(mol2)}")

# Isoquinoline core
iq_patt = Chem.MolFromSmarts("c1ccc2c(c1)ccnc2")
print(f"Isoquinoline match: {mol2.HasSubstructMatch(iq_patt)}")

if mol2.HasSubstructMatch(iq_patt):
    matches = mol2.GetSubstructMatches(iq_patt)
    for match in matches:
        print(f"  Match atoms: {match}")
        # Check what's attached outside the match
        for atom_idx in match:
            atom = mol2.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in match:
                    print(f"    Atom {atom_idx}({atom.GetSymbol()}) has outside neighbor: atom {nbr.GetIdx()}({nbr.GetSymbol()})")

# Actually run the retro reaction
rxn = AllChem.ReactionFromSmarts("c1ccc2c(c1)ccnc2>>NCCc1ccccc1.C=O")
ps = rxn.RunReactants((mol2,))
if ps and ps[0]:
    products = [Chem.MolToSmiles(p) for p in ps[0]]
    print(f"Retro products: {products}")