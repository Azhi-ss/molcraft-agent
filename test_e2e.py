#!/usr/bin/env python3
"""End-to-end test: SMILES -> 3D -> OpenBabel PDBQT -> Vina docking."""
import sys, os, tempfile
sys.path.insert(0, '/home/dministrator/Lab/clones/molcraft-agent/src')

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from openbabel import openbabel
from vina import Vina
from config import RECEPTOR_PDBQT, DOCKING_CENTER, DOCKING_SIZE

rdBase.DisableLog('rdApp.error')

smiles = 'Cc1nc(N)nc2n(C)c(C)c(C)c12'  # BMS-911543
print(f"Test SMILES: {smiles}")

# Step 1: SMILES -> RDKit Mol
mol = Chem.MolFromSmiles(smiles)
assert mol is not None, "MolFromSmiles failed"
print(f"Step 1 OK: {mol.GetNumAtoms()} heavy atoms")

# Step 2: Add H + 3D embed + optimize
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

# Step 3: 3D Mol -> PDBQT via OpenBabel (MolBlock string, no temp SDF file)
mol_block = Chem.MolToMolBlock(mol)
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("mol", "pdbqt")
# Add rigid receptor option (consistent with receptor.py)
obConversion.AddOption("r", openbabel.OBConversion.OUTOPTIONS)
obmol = openbabel.OBMol()
success = obConversion.ReadString(obmol, mol_block)
assert success, "OpenBabel ReadString failed"
print(f"Step 3 OK: OpenBabel read MolBlock, {obmol.NumAtoms()} atoms")

# Write PDBQT to temp file for Vina
fd, lig_pdbqt = tempfile.mkstemp(suffix=".pdbqt")
os.close(fd)
obConversion.WriteFile(obmol, lig_pdbqt)
with open(lig_pdbqt) as f:
    content = f.read()
    n_atoms = sum(1 for line in content.split('\n') if line.startswith('ATOM'))
print(f"Step 3b OK: PDBQT file written, {n_atoms} ATOM records")

# Step 4: Vina docking
v = Vina(sf_name="vina", seed=42, verbosity=0)
v.set_receptor(RECEPTOR_PDBQT)
v.set_ligand_from_file(lig_pdbqt)
v.compute_vina_maps(center=DOCKING_CENTER, box_size=DOCKING_SIZE)
v.dock(exhaustiveness=8, n_poses=5)
energies = v.energies(n_poses=1)
best_energy = float(energies[0][0])
print(f"Step 4 OK: Vina docking successful, binding energy = {best_energy:.3f} kcal/mol")

# Clean up
os.remove(lig_pdbqt)
print("\nAll steps passed! End-to-end pipeline works correctly.")