"""
Debug isoquinoline - use correct SMILES
"""
import sys
sys.path.insert(0, '/home/dministrator/Lab/clones/molcraft-agent/src')
from rdkit import Chem
from rdkit.Chem import AllChem

# Test various SMILES for isoquinoline
for smi in ["c1cncc2ccccc12", "c1ccc2ccnc2c1", "c1nc2ccccc2cc1"]:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        print(f"  {smi:25s} -> canonical: {Chem.MolToSmiles(mol):25s} atoms: {[a.GetSymbol() for a in mol.GetAtoms()]}")
    else:
        print(f"  {smi:25s} -> INVALID")

# Use the correct SMILES
target_smi = "c1cncc2ccccc12"  # isoquinoline, canonical
mol = Chem.MolFromSmiles(target_smi)
print(f"\nTarget: {Chem.MolToSmiles(mol)}")
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
    
    target_elems = {a.GetSymbol() for a in mol.GetAtoms()}
    reactant_elems = set()
    for r_smi in reactants:
        r_mol = Chem.MolFromSmiles(r_smi)
        if r_mol:
            reactant_elems |= {a.GetSymbol() for a in r_mol.GetAtoms()}
    print(f"Target elements: {target_elems}")
    print(f"Reactant elements: {reactant_elems}")
    suspicious = target_elems - reactant_elems - {'C', 'H', 'O', 'N'}
    print(f"Suspicious: {suspicious}")
    
    target_heavy = mol.GetNumHeavyAtoms()
    reactant_heavy = sum(Chem.MolFromSmiles(r).GetNumHeavyAtoms() for r in reactants if Chem.MolFromSmiles(r))
    ratio = target_heavy / reactant_heavy if reactant_heavy > 0 else 0
    print(f"Ratio: {target_heavy}/{reactant_heavy} = {ratio:.2f} (need 0.75-1.25)")

# Now test the full module
print("\n--- plan_synthesis_v2 on isoquinoline ---")
from synthesis_v2 import plan_synthesis_v2
result = plan_synthesis_v2(target_smi)
print(f"Result: {result}")

# Also test other key molecules
print("\n--- plan_synthesis_v2 on simple quinazoline ---")
result = plan_synthesis_v2("c1ccc2ncncc2c1")
print(f"Result: {result}")