#!/usr/bin/env python3
"""Test Meeko PDBQT after installing gemmi."""
import sys, os, tempfile
sys.path.insert(0, '/home/dministrator/Lab/clones/molcraft-agent/src')

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina
from config import RECEPTOR_PDBQT, DOCKING_CENTER, DOCKING_SIZE

rdBase.DisableLog('rdApp.error')

smiles = 'Cc1nc(N)nc2n(C)c(C)c(C)c12'  # BMS-911543
print(f"Test SMILES: {smiles}")

# Step 1: SMILES -> RDKit Mol -> 3D conformer
mol = Chem.MolFromSmiles(smiles)
assert mol is not None, "MolFromSmiles failed"
print(f"Step 1 OK: {mol.GetNumAtoms()} heavy atoms")

mol = Chem.AddHs(mol)
ret = AllChem.EmbedMolecule(mol, randomSeed=42)
if ret != 0:
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    ret = AllChem.EmbedMolecule(mol, params)
assert ret == 0, f"EmbedMolecule failed: {ret}"
print(f"Step 2 OK: {mol.GetNumAtoms()} atoms with H, 3D embedded")

try:
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
except Exception:
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
print("Step 2b: optimization done")

# Step 3: Meeko PDBQT conversion
preparator = MoleculePreparation()
setup_list = preparator.prepare(mol)
assert setup_list, "Meeko prepare returned empty list"
pdbqt_string = PDBQTWriterLegacy.write_string(setup_list[0])[0]
assert pdbqt_string and pdbqt_string.strip(), "Meeko PDBQT string is empty"
print(f"Step 3 OK: Meeko PDBQT, length: {len(pdbqt_string)}")
lines = pdbqt_string.split('\n')
for line in lines[:8]:
    print(f"  {line}")

# Write PDBQT to temp file
fd, lig_pdbqt = tempfile.mkstemp(suffix=".pdbqt")
os.close(fd)
with open(lig_pdbqt, "w") as f:
    f.write(pdbqt_string)
print(f"Step 3b OK: PDBQT file written to {lig_pdbqt}")

# Step 4: Vina docking
v = Vina(sf_name="vina", seed=42, verbosity=0)
v.set_receptor(RECEPTOR_PDBQT)
v.set_ligand_from_file(lig_pdbqt)
v.compute_vina_maps(center=DOCKING_CENTER, box_size=DOCKING_SIZE)
v.dock(exhaustiveness=8, n_poses=5)
energies = v.energies(n_poses=1)
best_energy = float(energies[0][0])
print(f"Step 4 OK: Vina docking successful, binding energy = {best_energy:.3f} kcal/mol")

os.remove(lig_pdbqt)
print("\nAll steps passed! Meeko pipeline works with gemmi installed.")