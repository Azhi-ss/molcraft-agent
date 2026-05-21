"""
Debug isoquinoline planning failure
"""
import sys
sys.path.insert(0, '/home/dministrator/Lab/clones/molcraft-agent/src')
from rdkit import Chem
from rdkit.Chem import AllChem

# Test: does the isoquinoline retro rule work on simple isoquinoline?
target_smi = "c1ccc2ccnc2c1"
mol = Chem.MolFromSmiles(target_smi)
print(f"Target: {target_smi}")
print(f"Atoms: {[a.GetSymbol() for a in mol.GetAtoms()]}")
print(f"Elements: {set(a.GetSymbol() for a in mol.GetAtoms())}")

# Test the isoquinoline rule
smarts_pat = "c1ccc2c(c1)ccnc2"
retro_smarts = "c1ccc2c(c1)ccnc2>>NCCc1ccccc1.C=O"

patt = Chem.MolFromSmarts(smarts_pat)
print(f"SMARTS match: {mol.HasSubstructMatch(patt)}")

rxn = AllChem.ReactionFromSmarts(retro_smarts)
ps = rxn.RunReactants((mol,))
if ps and ps[0]:
    reactants = [Chem.MolToSmiles(p) for p in ps[0]]
    print(f"Reactants: {reactants}")
    
    # Check element conservation
    target_elems = {a.GetSymbol() for a in mol.GetAtoms()}
    reactant_elems = set()
    for r_smi in reactants:
        r_mol = Chem.MolFromSmiles(r_smi)
        if r_mol:
            reactant_elems |= {a.GetSymbol() for a in r_mol.GetAtoms()}
    print(f"Target elements: {target_elems}")
    print(f"Reactant elements: {reactant_elems}")
    suspicious = target_elems - reactant_elems - {'C', 'H', 'O', 'N'}
    print(f"Suspicious (should be empty): {suspicious}")
    
    # Check ratio
    target_heavy = mol.GetNumHeavyAtoms()
    reactant_heavy = sum(Chem.MolFromSmiles(r).GetNumHeavyAtoms() for r in reactants if Chem.MolFromSmiles(r))
    ratio = target_heavy / reactant_heavy if reactant_heavy > 0 else 0
    print(f"Ratio: {target_heavy}/{reactant_heavy} = {ratio:.2f} (need 0.75-1.25)")
else:
    print("No reaction products")
    print(f"Rxn num products: {rxn.GetNumProductTemplates()}")
    print(f"Rxn num reactants: {rxn.GetNumReactantTemplates()}")

# Also test: does plan_synthesis work?
print("\n--- plan_synthesis_v2 on isoquinoline ---")
from synthesis_v2 import plan_synthesis_v2, plan_synthesis_recursive
result = plan_synthesis_recursive(target_smi, max_depth=1)
print(f"Result (max_depth=1): {result}")

result2 = plan_synthesis_v2(target_smi)
print(f"Result (v2, max_depth=3): {result2}")

# Test with explicit atoms
print("\n--- Multiple isoquinoline variants ---")
for smi in ["c1ccc2ccnc2c1", "c1cncc2ccccc12"]:
    mol = Chem.MolFromSmiles(smi)
    print(f"\n  {smi} -> {Chem.MolToSmiles(mol)}")
    print(f"  Atoms: {[a.GetSymbol() for a in mol.GetAtoms()]}")
    
    # Run through all rules
    for sp, rs in [
        ("c1ccc2c(c1)ccnc2", "c1ccc2c(c1)ccnc2>>NCCc1ccccc1.C=O"),
    ]:
        pat = Chem.MolFromSmarts(sp)
        if mol.HasSubstructMatch(pat):
            r = AllChem.ReactionFromSmarts(rs)
            p = r.RunReactants((mol,))
            if p and p[0]:
                prods = [Chem.MolToSmiles(pp) for pp in p[0]]
                print(f"  Matches rule '{sp}' -> {prods}")