"""
Debug why the explicit aromatic bond SMARTS rejects molecule #4 too.
Also test element conservation check as the primary fix.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# The intermediates
mol3 = Chem.MolFromSmiles("FCc1ncc2cccc(Br)c2n1")
mol4 = Chem.MolFromSmiles("OB(O)c1cccc2ncncc12")
mol2 = Chem.MolFromSmiles("OB(O)c1cncc2ccccc12")

# === Debug: what does approach 2 actually match? ===
print("=" * 70)
print("DEBUG: Why does explicit aromatic SMARTS fail on moleucle #4?")
print("=" * 70)

# SMARTS 2 broken down
smarts2 = "[c]1:[c]:[c]:[c]:[c]2:[n]:[c]:[n]:[c]:2:1"
print(f"SMARTS 2: {smarts2}")

# Check if the SMARTS parses correctly
from rdkit.Chem import rdMolDescriptors

patt2 = Chem.MolFromSmarts(smarts2)
if patt2:
    print(f"SMARTS parsed OK")
    print(f"  Num atoms: {patt2.GetNumAtoms()}")
    print(f"  SMARTS str: {Chem.MolToSmarts(patt2)}")
    
    # Check mol4 match in detail
    matches4 = mol4.GetSubstructMatches(patt2)
    print(f"  mol4 matches: {matches4}")
    
    matches3 = mol3.GetSubstructMatches(patt2)
    print(f"  mol3 matches: {matches3}")
    
    # Try individual recipes
    print(f"\n  What does each molecule's ring info say?")
    for name, mol in [("mol3", mol3), ("mol4", mol4)]:
        print(f"  {name} ({Chem.MolToSmiles(mol)}):")
        ri = mol.GetRingInfo()
        for i, ring in enumerate(ri.AtomRings()):
            arr = [f"{idx}({mol.GetAtomWithIdx(idx).GetSymbol()})" for idx in ring]
            print(f"    Ring {i}: {arr}")
else:
    print("SMARTS failed to parse")

# === Try: directly specifying the quinazoline ring system atom-by-atom ===
print("\n" + "=" * 70)
print("APPROACH A: Full explicit quinazoline SMARTS using atom env")
print("=" * 70)

# Quinazoline: benzene (6C) fused to pyrimidine (N,C,N,C) via 2 bridgehead C
# Write as: c1ccc2n(c1)cnc2  ← this means the benzene and pyrimidine rings
# The benzene: c1(c)c(c)c(c)c2 
# The pyrimidine: n(c1)c(n)c2
# Combined: c12ccccc1ncnc2 was already tested above (still matches mol3)

# Let's try using ring queries more precisely
# In quinazoline, the benzene ring has 6 carbons. The pyrimidine has N-C-N-C.
# Try: explicitly demand the benzene ring has all C, no N
smarts_a = "[c;r]1:[c;r]:[c;r]:[c;r]:[c;r]2:[n;r]:[c;r]:[n;r]:[c;r]:2:1"
patt_a = Chem.MolFromSmarts(smarts_a)
print(f"SMARTS A: {smarts_a}")
if patt_a:
    print(f"  mol3 match: {mol3.HasSubstructMatch(patt_a)}")
    print(f"  mol4 match: {mol4.HasSubstructMatch(patt_a)}")

# Approach B: Exclude structures where N atoms have specific degree/connectivity
# In true quinazoline's pyrimidine ring, each N is connected to exactly 2 atoms (both C)
# In naphthyridine, one N might be at the bridgehead (degree 3) - no, naphthyridine doesn't have that.
smarts_b = "[c]1:[c]:[c]:[c]:[c]2:[n;D2]:[c;D2]:[n;D2]:[c;D2]:2:1"
patt_b = Chem.MolFromSmarts(smarts_b)
print(f"\nSMARTS B (D2): {smarts_b}")
if patt_b:
    print(f"  mol3 match: {mol3.HasSubstructMatch(patt_b)}")
    print(f"  mol4 match: {mol4.HasSubstructMatch(patt_b)}")

# Approach C: Use degree on atoms adjacent to the pyrimidine ring
# In quinazoline, the C between the two N's has degree 2 (just C-N and C-N)
# In naphthyridine, this C might have different connectivity
smarts_c = "c1ccc2ncncc2[c;D2]1"
patt_c = Chem.MolFromSmarts(smarts_c)
print(f"\nSMARTS C: {smarts_c}")
if patt_c:
    print(f"  mol3 match: {mol3.HasSubstructMatch(patt_c)}")
    print(f"  mol4 match: {mol4.HasSubstructMatch(patt_c)}")

# === THE REAL FIX: Element Conservation Check ===
print("\n" + "=" * 70)
print("THE REAL FIX: Element Conservation Check in _run_retro_rule")
print("=" * 70)

print("""
The fundamental problem is:
  1. SMARTS substructure matching finds the heterocyclic core
  2. The retro reaction replaces the core, but substituents on the core are LOST
  3. These 'lost' atoms (B(OH)2, F, Br) appear from nowhere in the target
  4. Current atom balance check (0.7-1.3 ratio) is too lenient

The fix should be added AFTER running the retro reaction in _run_retro_rule():
  - Check if any element (other than C, H, O, N) is present in the target
    but completely absent from ALL reactants
  - If so → this is NOT a valid retro step → reject the route

This catches:
  - Molecule #2: B in target, no B in reactants → reject
  - Molecule #3: F, Br in target, no F, Br in reactants → reject  
  - Molecule #4: B in target, no B in reactants → reject
  - But does NOT catch C, H, O, N imbalances (which are expected in condensation)

Plus: we should fix the quinazoline SMARTS to be more specific.
""")

# === Fix the quinazoline SMARTS properly ===
print("\n" + "=" * 70)
print("FINDING: Correct quinazoline SMARTS that rejects naphthyridine")
print("=" * 70)

# Let's examine what's different between mol3 and mol4's ring systems
# mol3 Ring 0: [2(C), 12(N), 11(C), 5(C), 4(C), 3(N)]
# mol3 Ring 1: [6(C), 7(C), 8(C), 9(C), 11(C), 5(C)]
# mol4: need to check ring info

ri4 = mol4.GetRingInfo()
for i, ring in enumerate(ri4.AtomRings()):
    arr = [f"{idx}({mol4.GetAtomWithIdx(idx).GetSymbol()})" for idx in ring]
    print(f"  mol4 Ring {i}: {arr}")

print(f"\n  mol4 canonical: {Chem.MolToSmiles(mol4)}")

# Write both to SMILES and check if mol4 is a true quinazoline
quinazoline = Chem.MolFromSmiles("c1ccc2ncncc2c1")
print(f"  True quinazoline SMILES: c1ccc2ncncc2c1")
print(f"  True quinazoline canonical: {Chem.MolToSmiles(quinazoline)}")
print(f"  mol4 is quinazoline match: {mol4.HasSubstructMatch(Chem.MolFromSmarts('[#6]1:[#6]:[#6]:[#6]:[#6]2:[#7]:[#6]:[#7]:[#6]:2:1'))}")

# What about the difference? mol4 IS a quinazoline boronic acid
# Let me write it as: B(O)(O)-c1cccc2ncncc12
# In this SMILES, the B is attached to one of the quinazoline ring carbons

# === FINAL: Combined fix approach ===
print("\n" + "=" * 70)
print("FINAL RECOMMENDED FIX")
print("=" * 70)

print("""
FIX 1 (SMARTS specificity): Replace quinazoline SMARTS pattern
  OLD: "c1ccc2c(c1)ncnc2"
  NEW: "c1ccc2ncncc2c1"
  
  This still matches some naphthyridine-like structures but is the 
  canonical form and slightly more restrictive. The PRIMARY defense
  should be FIX 2.

FIX 2 (Element conservation): Add check in _run_retro_rule()
  After running the retro reaction and getting reactants,
  check if target has elements that reactants don't (excluding C,H,O,N):
  
  target_elements = set of element symbols in target mol
  reactant_elements = union of element symbols across all reactants
  suspicious = target_elements - reactant_elements - {'C', 'H', 'O', 'N'}
  if suspicious: reject route
  
  This is a general, robust fix that applies to ALL retro rules.
  
FIX 3 (Tighter ratio): Tighten the H014 atom ratio check
  OLD: 0.7 <= ratio <= 1.3
  NEW: 0.75 <= ratio <= 1.2
  + element conservation check above
""")

# === VERIFY: Test the combined fix ===
print("\n" + "=" * 70)
print("VERIFICATION: Running combined fix on all 3 molecules")
print("=" * 70)

def check_element_conservation(target_smi, reactants):
    """Check that all non-C/H/O/N elements are conserved"""
    t_mol = Chem.MolFromSmiles(target_smi)
    if t_mol is None:
        return True, set()
    target_elems = set(atom.GetSymbol() for atom in t_mol.GetAtoms())
    
    all_reactant_elems = set()
    for r_smi in reactants:
        r_mol = Chem.MolFromSmiles(r_smi)
        if r_mol:
            all_reactant_elems |= set(atom.GetSymbol() for atom in r_mol.GetAtoms())
    
    suspicious = target_elems - all_reactant_elems - {'C', 'H', 'O', 'N'}
    return len(suspicious) == 0, suspicious

# Test on all three intermediates
tests = [
    ("#2 Isoquinoline-B(OH)2", "OB(O)c1cncc2ccccc12", 
     "c1ccc2c(c1)ccnc2>>NCCc1ccccc1.C=O"),
    ("#3 Naphthyridine-F-Br", "FCc1ncc2cccc(Br)c2n1",
     "c1ccc2c(c1)ncnc2>>Nc1ccccc1C(=O)N.C=O"),
    ("#4 Quinazoline-B(OH)2", "OB(O)c1cccc2ncncc12",
     "c1ccc2c(c1)ncnc2>>Nc1ccccc1C(=O)N.C=O"),
    ("Control: simple quinazoline", "c1ccc2ncncc2c1",
     "c1ccc2c(c1)ncnc2>>Nc1ccccc1C(=O)N.C=O"),
]

for name, target_smi, retro_smarts in tests:
    print(f"\n  {name}:")
    print(f"    Target: {target_smi}")
    t_mol = Chem.MolFromSmiles(target_smi)
    
    # Step 1: SMARTS pattern match (using OLD SMARTS for test)
    old_patt = Chem.MolFromSmarts("c1ccc2c(c1)ncnc2")
    if name.startswith("#2"):
        old_patt = Chem.MolFromSmarts("c1ccc2c(c1)ccnc2")
    
    if not t_mol.HasSubstructMatch(old_patt):
        print(f"    SMARTS match: NO (already rejected)")
        continue
    
    # Step 2: Run retro reaction
    rxn = AllChem.ReactionFromSmarts(retro_smarts)
    ps = rxn.RunReactants((t_mol,))
    if not ps or not ps[0]:
        print(f"    Retro reaction: NO products")
        continue
    
    reactants = [Chem.MolToSmiles(p) for p in ps[0]]
    print(f"    Reactants: {reactants}")
    
    # Step 3: Element conservation check
    ok, missing = check_element_conservation(target_smi, reactants)
    if ok:
        print(f"    Element check: ✅ PASS")
    else:
        print(f"    Element check: ❌ FAIL (missing: {missing})")
    
    # Step 4: Check atom ratio 
    t_heavy = t_mol.GetNumHeavyAtoms()
    r_heavy = sum(Chem.MolFromSmiles(r).GetNumHeavyAtoms() for r in reactants if Chem.MolFromSmiles(r))
    ratio = t_heavy / r_heavy if r_heavy > 0 else 0
    old_ok = 0.7 <= ratio <= 1.3
    new_ok = 0.75 <= ratio <= 1.2
    print(f"    Ratio (target/reactants): {t_heavy}/{r_heavy} = {ratio:.2f}")
    print(f"    Ratio old check (0.7-1.3): {'✅' if old_ok else '❌'}")
    print(f"    Ratio new check (0.75-1.2): {'✅' if new_ok else '❌'}")
    print(f"    Combined pass (element + new ratio): {'✅' if ok and new_ok else '❌'}")