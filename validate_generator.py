"""验证 generator.py 中所有新增 SMILES/SMARTS 的 RDKit 可解析性。"""
from rdkit import Chem

def validate_smiles_list(name, items):
    failed = []
    for s in items:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            failed.append(s)
    print(f"{name}: {len(items)} total, {len(failed)} failed")
    if failed:
        for f in failed:
            print(f"  FAILED: {f}")
    return failed

def validate_smarts_list(name, items):
    failed = []
    for s in items:
        mol = Chem.MolFromSmarts(s)
        if mol is None:
            failed.append(s)
    print(f"{name}: {len(items)} total, {len(failed)} failed")
    if failed:
        for f in failed:
            print(f"  FAILED: {f}")
    return failed

# 取代基
substituents = [
    'F', 'Cl', '[OH]', '[NH2]', 'C',
    'C#N', 'CO', 'COC', 'CON',
    'OCCN', 'NC(=O)C', 'NC(=O)',
    'C(F)(F)F', 'C(F)(F)F.C', 'c1cccc(F)c1', 'c1cccc(Cl)c1',
    'C1CNC1', 'C1CCNC1', 'C1CCNCC1', 'C1CNCCN1', 'C1COCCN1',
    'C1CC1', 'C1CCC1',
    'c1ccncc1', 'c1cncnc1', 'c1cnc[nH]1', 'c1ccnc1',
    '[N+](=O)[O-]', 'S(=O)(=O)C', 'S', 'C(=O)OC', 'C(=O)O',
]
f1 = validate_smiles_list("Substituents", substituents)

# 新增骨架
new_scaffolds = [
    'C1CNCCN1C', 'C1CN(C)CCN1',
    'c1cc[nH]n1', 'c1c[nH]ncc1C', 'c1cnc[nH]c1',
    'c1ccc2c(c1)[nH]n2', 'c1ccc2c(c1)n[nH]2',
    'c1cnc2[nH]cc[nH]2c1', 'c1cnc2[nH]ccc2c1', 'c1cnc2cc[nH]2c1',
    'c1ccc2ncnc2c1', 'c1ccc2c(c1)nccn2',
    'c1cnc2ncncc2n1', 'c1cnc2cncnc2n1',
    'c1cscn1', 'c1cscnc1', 'c1nccs1', 'c1ncncs1', 'c1nccns1',
    'c1cncn1', 'c1nnnc1', 'c1nncn1',
    'c1ccc2c(c1)nnc2', 'c1ccc2c(c1)ncns2',
    'c1ccc2c(c1)CCNC2', 'c1ccc2c(c1)CCN(C)C2', 'c1ccc2c(c1)CCCC2',
    'c1ccc2c(c1)CCCN2', 'c1ccc2ocnc2c1',
    'c1cnccn1', 'c1ccnnc1', 'c1cnccnc1N',
]
f2 = validate_smiles_list("New scaffolds", new_scaffolds)

# SMARTS
smarts_list = ['[n]', '[nH]', '[nH0]', 'n1ccnc1', '[nH]c2ncnc2']
f3 = validate_smarts_list("SMARTS", smarts_list)

# 激酶优先取代基
kinase_priority = [
    'C#N', '[NH2]', '[OH]', 'c1ccncc1', 'c1cncnc1',
    'c1cnc[nH]1', 'c1cc[nH]n1', 'F', 'Cl',
    'C1CC1', 'C(F)(F)F', 'CO', 'NC(=O)',
    'OCCN', 'C1CNCCN1', 'C1COCCN1',
]
f4 = validate_smiles_list("Kinase priority subs", kinase_priority)

# 验证全部 SCAFFOLDS（包括原始的）
original_scaffolds = [
    "c1ccc(cc1)", "c1ccccc1C", "c1ccc(cc1)O", "c1ccc(cc1)N", "COc1ccc(cc1)", "c1ccc(cc1)CN",
    "c1ccc2ccccc2c1", "c1ccc2c(c1)cccc2",
    "c1ccc(nc1)", "c1cncnc1", "c1cnccc1", "c1c[nH]cn1",
    "C1CCOC1", "C1CCNC1", "C1CCNCC1", "C1CNCCN1", "C1COCCN1", "C1CSCN1", "C1COC1", "C1CNC1",
    "c1ccc2cncnc2c1", "c1ccc2nc[nH]c2c1", "c1ccc2c(c1)cccn2", "c1ccc2c(c1)ocn2",
    "c1ccc2c(c1)scn2", "c1ccc2ncccc2c1", "c1ccc2ccncc2c1", "c1ccc2CCNc2c1",
    "c1ccc2c(c1)CCN2", "c1ccc2c(c1)CCNC2", "c1ccc2c(c1)N=CN2", "c1nc2c([nH]1)cccc2",
    "c1[nH]cnc2ncnc12", "c1cnc2ncncc2n1", "c1ccc2c(c1)cncn2", "c1ccc2c(c1)ncnc2",
    "c1ccc2c(c1)OCCO2", "c1ccc2c(c1)CCO2",
    "c1ccc(cc1)C(=O)O", "c1ccc(cc1)C(=O)N", "c1ccc(cc1)S(=O)(=O)N",
    "c1ccc(cc1)C(=O)Nc2ccccc2", "c1ccc(cc1)NC(=O)c2ccccc2", "Cc1ccc(cc1)S(=O)(=O)Nc2ccccc2",
    "c1cc(ccc1F)F", "c1cc(ccc1Cl)Cl", "c1ccc(cc1)CC(=O)O",
    "c1ccc(cc1)C(=O)c2ccccc2", "c1ccc(cc1)OCc2ccccc2",
    "c1ccc2c(c1)CCCC2", "c1ccc2c(c1)CCCCC2", "c1ccc2c(c1)NCCC2",
    "c1ccc2c(c1)OCC2", "c1ccc2c(c1)SCC2", "c1ccc2c(c1)CCO2",
    "c1nc2c(n1)ncn2", "c1cnc2[nH]ccc2c1", "c1cc2nccc2n1", "c1cc2ncnc2n1", "c1ccc2cccn2c1",
    "C1CC2CCC1C2", "C1CN2CCC1CC2", "C1CC2CCC(C1)N2", "C1C2CC1C2",
    "C1CC2(CCNCC2)NC1", "c1ccc2c(c1)CC3(CCCCC3)N2", "C1NCC2(COC2)C1",
    "C1CCCNCC1", "C1CCCOCC1", "O=S1(=O)CCNCC1",
]
f5 = validate_smiles_list("Original scaffolds", original_scaffolds)

# 验证铰链骨架
hinge_scaffolds = [
    "c1cnc2[nH]ccc2c1", "c1[nH]cnc2ncnc12", "c1cc2ncnc2n1", "c1cc2nccc2n1",
    "c1ccc2c(c1)ncnc2", "c1ccc2c(c1)cccn2", "c1cncnc1", "c1nc2c(n1)ncn2",
    "c1ccc2c(c1)[nH]n2", "c1ccc2c(c1)n[nH]2",
    "c1cc[nH]n1", "c1cnc2[nH]cc[nH]2c1", "c1cnc2cc[nH]2c1",
    "c1ccc2c(c1)nccn2", "c1cnc2ncncc2n1", "c1cnc2cncnc2n1",
    "c1cscn1", "c1cncn1",
    "c1ccc2c(c1)CCNC2", "c1ccc2c(c1)CCN(C)C2", "c1cnccnc1N",
]
f6 = validate_smiles_list("Hinge scaffolds", hinge_scaffolds)

# 汇总
all_failed = f1 + f2 + f3 + f4 + f5 + f6
if all_failed:
    print(f"\nTOTAL FAILED: {len(all_failed)}")
    for f in all_failed:
        print(f"  {f}")
else:
    print(f"\nALL PASSED - 0 failures")

# 尝试导入 generator
try:
    import sys
    sys.path.insert(0, 'src')
    from generator import SCAFFOLDS, KINASE_HINGE_SCAFFOLDS, _add_substituent, _mutate_mol
    print(f"\nImport successful!")
    print(f"  SCAFFOLDS: {len(SCAFFOLDS)} entries")
    print(f"  KINASE_HINGE_SCAFFOLDS: {len(KINASE_HINGE_SCAFFOLDS)} entries")
    
    # 验证所有 SCAFFOLDS 可解析
    scaffold_failed = []
    for s in SCAFFOLDS:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            scaffold_failed.append(s)
    print(f"  SCAFFOLDS parse check: {len(SCAFFOLDS)} total, {len(scaffold_failed)} failed")
    
    hinge_failed = []
    for s in KINASE_HINGE_SCAFFOLDS:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            hinge_failed.append(s)
    print(f"  HINGE parse check: {len(KINASE_HINGE_SCAFFOLDS)} total, {len(hinge_failed)} failed")
    
    # 验证新增函数
    mol = Chem.MolFromSmiles('c1ccncc1')
    result = _add_substituent(mol)
    if result is not None:
        print(f"  _add_substituent test: OK (produced {Chem.MolToSmiles(result)})")
    else:
        print(f"  _add_substituent test: returned None (may be ok for some cases)")
    
    # 验证位置感知变异
    mol2 = Chem.MolFromSmiles('c1cnc2[nH]ccc2c1')  # 7-azaindole
    result2 = _mutate_mol(mol2, docking_guidance={"priority_weight": 0.8})
    if result2 is not None:
        print(f"  _mutate_mol_position_aware test: OK (produced {Chem.MolToSmiles(result2)})")
    else:
        print(f"  _mutate_mol_position_aware test: returned None")
    
    # 验证无 docking_guidance 的默认行为
    result3 = _mutate_mol(mol)
    if result3 is not None:
        print(f"  _mutate_mol (no guidance) test: OK")
    else:
        print(f"  _mutate_mol (no guidance) test: returned None")
    
except Exception as e:
    print(f"\nImport FAILED: {e}")