"""
Deep analysis of SMARTS matching and fixes for molecules #2, #3, #4
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# === STEP 1: Understand WHY quinazoline SMARTS matches non-quinazolines ===

print("=" * 80)
print("STEP 1: Why does quinazoline SMARTS match molecule #3?")
print("=" * 80)

# The problem molecules
mol3 = Chem.MolFromSmiles("FCc1ncc2cccc(Br)c2n1")
mol4 = Chem.MolFromSmiles("OB(O)c1cccc2ncncc12")

# Current quinazoline SMARTS
old_qz = "c1ccc2c(c1)ncnc2"

print(f"\nMolecule #3: FCc1ncc2cccc(Br)c2n1")
print(f"  Canonical: {Chem.MolToSmiles(mol3)}")
print(f"  Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol3)}")
print(f"  Rings: {mol3.GetRingInfo().NumRings()}")
print(f"  Ring atoms: {[a.GetIdx() for a in mol3.GetAtoms() if a.IsInRing()]}")

# What IS molecule #3's ring system?
# Describe each ring
ri = mol3.GetRingInfo()
for i, ring in enumerate(ri.AtomRings()):
    atoms_str = ", ".join([f"{idx}({mol3.GetAtomWithIdx(idx).GetSymbol()})" for idx in ring])
    print(f"  Ring {i}: [{atoms_str}]")

print(f"\nMolecule #4: OB(O)c1cccc2ncncc12")
print(f"  Canonical: {Chem.MolToSmiles(mol4)}")
print(f"  Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol4)}")

# Test the current SMARTS on both
old_patt = Chem.MolFromSmarts(old_qz)
print(f"\nOld SMARTS '{old_qz}':")
print(f"  Molecule #3 (naphthyridine) match: {mol3.HasSubstructMatch(old_patt)}")
print(f"  Molecule #4 (quinazoline) match:   {mol4.HasSubstructMatch(old_patt)}")

# The match details
matches3 = mol3.GetSubstructMatches(old_patt)
for m in matches3:
    print(f"    Match atoms (mol #3): {m}")
    symbols = [f"{idx}({mol3.GetAtomWithIdx(idx).GetSymbol()})" for idx in m]
    print(f"    Symbols: {symbols}")

matches4 = mol4.GetSubstructMatches(old_patt)
for m in matches4:
    print(f"    Match atoms (mol #4): {m}")
    symbols = [f"{idx}({mol4.GetAtomWithIdx(idx).GetSymbol()})" for idx in m]
    print(f"    Symbols: {symbols}")

# === STEP 2: Test improved SMARTS patterns ===

print("\n" + "=" * 80)
print("STEP 2: Test improved SMARTS patterns for quinazoline")
print("=" * 80)

# Approach 1: Use the canonical SMILES form of quinazoline as SMARTS
# c1ccc2ncncc2c1 - this is the canonical form of quinazoline
improved_qz_1 = "c1ccc2ncncc2c1"
print(f"\nImproved SMARTS 1: '{improved_qz_1}' (canonical quinazoline)")
patt1 = Chem.MolFromSmarts(improved_qz_1)
print(f"  Molecule #3 match: {mol3.HasSubstructMatch(patt1)}")
print(f"  Molecule #4 match: {mol4.HasSubstructMatch(patt1)}")
if mol4.HasSubstructMatch(patt1):
    m4m = mol4.GetSubstructMatches(patt1)
    print(f"    Match atoms (mol #4): {m4m}")

# Approach 2: Use isotope labeling to force exact N positions
# In quinazoline, ring 2 has: N at position 1, C at 2, N at 3, C at 4 (from bridgehead)
# Using c:n:c:n: pattern explicitly with ring closures
improved_qz_2 = "[c]1:[c]:[c]:[c]:[c]2:[n]:[c]:[n]:[c]:2:1"
print(f"\nImproved SMARTS 2: '{improved_qz_2}' (explicit arom bonds)")
patt2 = Chem.MolFromSmarts(improved_qz_2)
if patt2:
    print(f"  Molecule #3 match: {mol3.HasSubstructMatch(patt2)}")
    print(f"  Molecule #4 match: {mol4.HasSubstructMatch(patt2)}")
else:
    print("  ERROR: Invalid SMARTS")

# Approach 3: Use Degree specification on bridgehead carbons
# In quinazoline, the bridgehead carbons connecting to the benzene side have degree 3
# and are specifically connected to N
improved_qz_3 = "c12ccccc1ncnc2"
print(f"\nImproved SMARTS 3: '{improved_qz_3}'")
patt3 = Chem.MolFromSmarts(improved_qz_3)
print(f"  Molecule #3 match: {mol3.HasSubstructMatch(patt3)}")
print(f"  Molecule #4 match: {mol4.HasSubstructMatch(patt3)}")
if mol4.HasSubstructMatch(patt3):
    print(f"    Match atoms (mol #4): {mol4.GetSubstructMatches(patt3)}")
if mol3.HasSubstructMatch(patt3):
    print(f"    Match atoms (mol #3): {mol3.GetSubstructMatches(patt3)}")

# Approach 4: Enforce degree on the N atoms - in quinazoline's pyrimidine ring
# the N atoms are degree 2 (connected to 2 carbons), not degree 1 (edge of ring)
improved_qz_4 = "c1ccc2[n;D2]c[n;D2]c2c1"
print(f"\nImproved SMARTS 4: '{improved_qz_4}' (D2 on N)")
patt4 = Chem.MolFromSmarts(improved_qz_4)
print(f"  Molecule #3 match: {mol3.HasSubstructMatch(patt4)}")
print(f"  Molecule #4 match: {mol4.HasSubstructMatch(patt4)}")

# Approach 5: Exclude ipso substitution - use !$(c(:[n])...) patterns
# This prevents matching when there's something on the other side
improved_qz_5 = "c1ccc2c(c1)[n;r]c[n;r]c2"
print(f"\nImproved SMARTS 5: '{improved_qz_5}'")
patt5 = Chem.MolFromSmarts(improved_qz_5)
print(f"  Molecule #3 match: {mol3.HasSubstructMatch(patt5)}")
print(f"  Molecule #4 match: {mol4.HasSubstructMatch(patt5)}")

# === STEP 3: Also test isoquinoline fix ===

print("\n" + "=" * 80)
print("STEP 3: Test improved isoquinoline SMARTS")
print("=" * 80)

mol2 = Chem.MolFromSmiles("OB(O)c1cncc2ccccc12")
print(f"\nMolecule #2 intermediate: OB(O)c1cncc2ccccc12")
print(f"  Canonical: {Chem.MolToSmiles(mol2)}")

old_iq = "c1ccc2c(c1)ccnc2"
improved_iq = "c1ccc2ccnc2c1"  # canonical form

print(f"Old SMARTS '{old_iq}': {mol2.HasSubstructMatch(Chem.MolFromSmarts(old_iq))}")
print(f"Improved SMARTS '{improved_iq}': {mol2.HasSubstructMatch(Chem.MolFromSmarts(improved_iq))}")


# === STEP 4: The ATOM BALANCE problem ===

print("\n" + "=" * 80)
print("STEP 4: Element consistency check - do reactants have all elements?")
print("=" * 80)

# For molecule #2's first step
target_smi = "OB(O)c1cncc2ccccc12"
target_mol = Chem.MolFromSmiles(target_smi)
rxn2 = AllChem.ReactionFromSmarts("c1ccc2c(c1)ccnc2>>NCCc1ccccc1.C=O")
ps2 = rxn2.RunReactants((target_mol,))
reactants2 = [Chem.MolToSmiles(p) for p in ps2[0]]

def get_elements(mol):
    """Get set of element symbols in a molecule"""
    return set(atom.GetSymbol() for atom in mol.GetAtoms())

target_elems = get_elements(target_mol)
all_reactant_elems = set()
for r_smi in reactants2:
    r_mol = Chem.MolFromSmiles(r_smi)
    if r_mol:
        all_reactant_elems |= get_elements(r_mol)

print(f"\nMolecule #2 first step:")
print(f"  Target: {target_smi}, elements: {target_elems}")
print(f"  Reactants: {reactants2}, elements: {all_reactant_elems}")
missing = target_elems - all_reactant_elems
print(f"  Elements present in target but absent in reactants: {missing}")

# For molecule #3's first step
target_smi3 = "FCc1ncc2cccc(Br)c2n1"
target_mol3 = Chem.MolFromSmiles(target_smi3)
rxn3 = AllChem.ReactionFromSmarts("c1ccc2c(c1)ncnc2>>Nc1ccccc1C(=O)N.C=O")
ps3 = rxn3.RunReactants((target_mol3,))
reactants3 = [Chem.MolToSmiles(p) for p in ps3[0]]

target_elems3 = get_elements(target_mol3)
all_reactant_elems3 = set()
for r_smi in reactants3:
    r_mol = Chem.MolFromSmiles(r_smi)
    if r_mol:
        all_reactant_elems3 |= get_elements(r_mol)

print(f"\nMolecule #3 first step:")
print(f"  Target: {target_smi3}, elements: {target_elems3}")
print(f"  Reactants: {reactants3}, elements: {all_reactant_elems3}")
missing3 = target_elems3 - all_reactant_elems3
print(f"  Elements present in target but absent in reactants: {missing3}")

# For molecule #4's first step
target_smi4 = "OB(O)c1cccc2ncncc12"
target_mol4 = Chem.MolFromSmiles(target_smi4)
rxn4 = AllChem.ReactionFromSmarts("c1ccc2c(c1)ncnc2>>Nc1ccccc1C(=O)N.C=O")
ps4 = rxn4.RunReactants((target_mol4,))
reactants4 = [Chem.MolToSmiles(p) for p in ps4[0]]

target_elems4 = get_elements(target_mol4)
all_reactant_elems4 = set()
for r_smi in reactants4:
    r_mol = Chem.MolFromSmiles(r_smi)
    if r_mol:
        all_reactant_elems4 |= get_elements(r_mol)

print(f"\nMolecule #4 first step:")
print(f"  Target: {target_smi4}, elements: {target_elems4}")
print(f"  Reactants: {reactants4}, elements: {all_reactant_elems4}")
missing4 = target_elems4 - all_reactant_elems4
print(f"  Elements present in target but absent in reactants: {missing4}")

# === STEP 5: Test the fix - retro reaction with element check ===

print("\n" + "=" * 80)
print("STEP 5: Run retro rules with BOTH SMARTS fix AND element check")
print("=" * 80)

def check_element_conservation(target_smi, reactants):
    """Check that all elements in target appear somewhere in reactants"""
    t_mol = Chem.MolFromSmiles(target_smi)
    if t_mol is None:
        return True
    target_elements = get_elements(t_mol)
    
    reactant_elements = set()
    for r_smi in reactants:
        r_mol = Chem.MolFromSmiles(r_smi)
        if r_mol:
            reactant_elements |= get_elements(r_mol)
    
    missing = target_elements - reactant_elements
    # C, H, O are expected to be lost in condensation reactions
    # But B, F, Cl, Br, I, S, P, Si should NOT appear from nowhere
    suspicious_missing = missing - {'C', 'H', 'O', 'N'}
    if suspicious_missing:
        return False, suspicious_missing
    return True, missing

# Test the improved quinazoline SMARTS
print("\nTesting improved quinazoline SMARTS (approach 1):")
improved_qz_smarts = "c1ccc2ncncc2c1"
improved_rxn = "c1ccc2ncncc2c1>>Nc1ccccc1C(=O)N.C=O"

# Molecule #3 with improved SMARTS
print(f"\n  Molecule #3 (naphthyridine) with improved SMARTS:")
patt = Chem.MolFromSmarts(improved_qz_smarts)
match3 = mol3.HasSubstructMatch(patt)
print(f"    SMARTS match: {match3}")

# Molecule #4 with improved SMARTS
print(f"\n  Molecule #4 (true quinazoline) with improved SMARTS:")
match4 = mol4.HasSubstructMatch(patt)
print(f"    SMARTS match: {match4}")

# Also check: does molecule #2 match quinazoline SMARTS?
print(f"\n  Molecule #2 (isoquinoline) with quinazoline SMARTS:")
match2_qz = mol2.HasSubstructMatch(patt)
print(f"    SMARTS match: {match2_qz}")