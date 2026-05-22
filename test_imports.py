#!/usr/bin/env python3
"""Test script: verify Meeko 0.7.1 works with RDKit 2026.3.1."""
import sys
import os
import tempfile
print(f"Python: {sys.version}")

try:
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem
    print(f"RDKit OK")
except Exception as e:
    print(f"RDKit error: {e}")

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    print(f"Meeko OK")
except Exception as e:
    print(f"Meeko error: {e}")

try:
    from openbabel import openbabel
    print("OpenBabel OK")
except Exception as e:
    print(f"OpenBabel error: {e}")

# Test full SMILES -> 3D -> PDBQT pipeline with OpenBabel
smiles = 'Cc1nc(N)nc2n(C)c(C)c(C)c12'  # BMS-911543
print(f"\nTest SMILES: {smiles}")

try:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("MolFromSmiles failed!")
    else:
        print(f"MolFromSmiles OK, {mol.GetNumAtoms()} heavy atoms")
        mol = Chem.AddHs(mol)
        print(f"AddHs OK, {mol.GetNumAtoms()} total atoms")
        
        ret = AllChem.EmbedMolecule(mol, randomSeed=42)
        print(f"EmbedMolecule ret: {ret}")
        if ret == 0:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            print("MMFF optimization OK")
            
            # Method 1: Try Meeko (should work now with 0.7.1)
            try:
                preparator = MoleculePreparation()
                setup_list = preparator.prepare(mol)
                pdbqt_string = PDBQTWriterLegacy.write_string(setup_list[0])[0]
                print(f"Meeko PDBQT OK, length: {len(pdbqt_string)}")
                print(f"  First line: {pdbqt_string.split(chr(10))[0]}")
            except Exception as e:
                print(f"Meeko PDBQT failed: {e}")
            
            # Method 2: OpenBabel via SDF
            try:
                # Write SDF
                fd, sdf_path = tempfile.mkstemp(suffix=".sdf")
                os.close(fd)
                writer = Chem.SDWriter(sdf_path)
                writer.write(mol)
                writer.close()
                
                # OpenBabel SDF -> PDBQT (no H added, no charges - consistent with receptor.py)
                obConversion = openbabel.OBConversion()
                obConversion.SetInAndOutFormats("sdf", "pdbqt")
                obConversion.AddOption("r", openbabel.OBConversion.OUTOPTIONS)  # rigid receptor
                obmol = openbabel.OBMol()
                obConversion.ReadFile(obmol, sdf_path)
                print(f"OpenBabel read SDF: {obmol.NumAtoms()} atoms")
                
                fd2, pdbqt_path = tempfile.mkstemp(suffix=".pdbqt")
                os.close(fd2)
                obConversion.WriteFile(obmol, pdbqt_path)
                
                with open(pdbqt_path) as f:
                    pdbqt_content = f.read()
                    print(f"OpenBabel PDBQT OK, length: {len(pdbqt_content)}")
                    lines = pdbqt_content.split('\n')
                    for line in lines[:5]:
                        print(f"  {line}")
                
                os.remove(sdf_path)
                os.remove(pdbqt_path)
            except Exception as e:
                print(f"OpenBabel SDF->PDBQT failed: {e}")
                import traceback
                traceback.print_exc()

            # Method 3: OpenBabel via MolBlock string (no temp file for SDF)
            try:
                mol_block = Chem.MolToMolBlock(mol)
                obConversion = openbabel.OBConversion()
                obConversion.SetInAndOutFormats("mol", "pdbqt")
                obConversion.AddOption("r", openbabel.OBConversion.OUTOPTIONS)
                obmol = openbabel.OBMol()
                obConversion.ReadString(obmol, mol_block)
                print(f"OpenBabel read MolBlock: {obmol.NumAtoms()} atoms")
                
                pdbqt_string_ob = obConversion.WriteString(obmol)
                print(f"OpenBabel MolBlock->PDBQT OK, length: {len(pdbqt_string_ob)}")
                lines = pdbqt_string_ob.split('\n')
                for line in lines[:5]:
                    print(f"  {line}")
            except Exception as e:
                print(f"OpenBabel MolBlock->PDBQT failed: {e}")
                import traceback
                traceback.print_exc()

except Exception as e:
    print(f"Pipeline error: {e}")
    import traceback
    traceback.print_exc()

# Test with Vina docking
print("\n--- Test Vina docking with OpenBabel PDBQT ---")
try:
    from vina import Vina
    from config import RECEPTOR_PDBQT, DOCKING_CENTER, DOCKING_SIZE
    print(f"Config: receptor={RECEPTOR_PDBQT}, center={DOCKING_CENTER}, size={DOCKING_SIZE}")
    print(f"Receptor exists: {os.path.exists(RECEPTOR_PDBQT)}")
    
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    
    # Generate PDBQT via OpenBabel
    mol_block = Chem.MolToMolBlock(mol)
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "pdbqt")
    obConversion.AddOption("r", openbabel.OBConversion.OUTOPTIONS)
    obmol = openbabel.OBMol()
    obConversion.ReadString(obmol, mol_block)
    
    fd, lig_pdbqt = tempfile.mkstemp(suffix=".pdbqt")
    os.close(fd)
    obConversion.WriteFile(obmol, lig_pdbqt)
    
    v = Vina(sf_name="vina", seed=42, verbosity=0)
    v.set_receptor(RECEPTOR_PDBQT)
    v.set_ligand_from_file(lig_pdbqt)
    v.compute_vina_maps(center=DOCKING_CENTER, box_size=DOCKING_SIZE)
    v.dock(exhaustiveness=8, n_poses=5)
    energies = v.energies(n_poses=1)
    best_energy = float(energies[0][0])
    print(f"Vina docking OK! Binding energy: {best_energy:.3f} kcal/mol")
    
    os.remove(lig_pdbqt)
except Exception as e:
    print(f"Vina docking failed: {e}")
    import traceback
    traceback.print_exc()