"""分子对接模块：使用 AutoDock Vina 计算结合自由能。"""
import os
import sys
import tempfile
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina
from config import RECEPTOR_PDBQT, DOCKING_CENTER, DOCKING_SIZE, DOCKING_EXHAUSTIVENESS
from receptor import prepare_receptor

# 压制 RDKit 错误输出
rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')


def smiles_to_pdbqt(smiles: str, output_path: str = None):
    """使用 RDKit + Meeko 将 SMILES 转换为 PDBQT 文件。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    ret = AllChem.EmbedMolecule(mol, randomSeed=42)
    if ret != 0:
        ret = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if ret != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass

    preparator = MoleculePreparation()
    setup_list = preparator.prepare(mol)
    if not setup_list:
        return None

    pdbqt_string = PDBQTWriterLegacy.write_string(setup_list[0])[0]

    if output_path is None:
        fd, output_path = tempfile.mkstemp(suffix=".pdbqt")
        os.close(fd)

    with open(output_path, "w") as f:
        f.write(pdbqt_string)

    return output_path


def dock_molecule(smiles: str, center=None, size=None, exhaustiveness=None, seed: int = 0):
    """对单个分子进行对接，返回结合能（kcal/mol）。

    参数:
        smiles: 配体 SMILES
        center, size: 对接盒子参数
        exhaustiveness: Vina 搜索精度
        seed: 随机种子（0 = 使用默认 42）

    返回字典包含:
        - binding_energy: 最佳构象能量（负值=越好）
        - poses: 构象能量列表
        - success: 是否成功
        - error: 失败原因
    """
    receptor = prepare_receptor()
    ligand_pdbqt = smiles_to_pdbqt(smiles)
    if ligand_pdbqt is None:
        return {"success": False, "error": "配体准备失败"}

    if center is None:
        center = DOCKING_CENTER
    if size is None:
        size = DOCKING_SIZE
    if exhaustiveness is None:
        exhaustiveness = DOCKING_EXHAUSTIVENESS

    vina_seed = seed if seed != 0 else 42
    try:
        v = Vina(sf_name="vina", seed=vina_seed, verbosity=0)
        v.set_receptor(receptor)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=center, box_size=size)
        v.dock(exhaustiveness=exhaustiveness, n_poses=5)
        energies = v.energies(n_poses=1)
        best_energy = float(energies[0][0])

        # 清理临时文件
        if ligand_pdbqt.startswith(tempfile.gettempdir()):
            os.remove(ligand_pdbqt)

        return {
            "success": True,
            "binding_energy": round(best_energy, 3),
            "poses": energies.tolist(),
        }
    except Exception as e:
        if ligand_pdbqt and ligand_pdbqt.startswith(tempfile.gettempdir()):
            try:
                os.remove(ligand_pdbqt)
            except Exception:
                pass
        return {"success": False, "error": str(e)}


def dock_molecule_consensus(smiles: str, center=None, size=None, exhaustiveness=None,
                            n_runs: int = 3, seeds: list = None):
    """共识对接：多次独立对接取中位数结合能（H009）。

    Coscientist 核心模式："performing experiments multiple times" —
    通过不同随机种子执行 N 次独立对接，消除单次运行的随机偏差，
    取中位数作为共识评分。

    返回:
        - success: 是否有足够成功的对接（>=2/3）
        - binding_energy: 中位数结合能（kcal/mol）
        - all_energies: 所有成功对接的能量列表
        - std_energy: 多次运行的能量标准差
        - n_success: 成功对接次数
    """
    import statistics
    if seeds is None:
        seeds = [42, 123, 456]
    seeds = seeds[:n_runs]

    energies = []
    success_count = 0
    for seed in seeds:
        result = dock_molecule(smiles, center, size, exhaustiveness, seed=seed)
        if result.get("success"):
            energies.append(result["binding_energy"])
            success_count += 1

    if success_count < 2:
        # 不足两次成功：降级为单次结果或失败
        if success_count == 1:
            return {
                "success": True,
                "binding_energy": energies[0],
                "all_energies": energies,
                "std_energy": 0.0,
                "n_success": 1,
            }
        return {"success": False, "error": f"共识对接失败: {success_count}/{n_runs}"}

    med_energy = round(statistics.median(energies), 3)
    std_energy = round(statistics.stdev(energies) if success_count > 1 else 0.0, 3)
    return {
        "success": True,
        "binding_energy": med_energy,
        "all_energies": energies,
        "std_energy": std_energy,
        "n_success": success_count,
    }


def batch_dock(molecules, center=None, size=None):
    """批量对接分子并返回结果。"""
    results = []
    total_count = len(molecules)
    completed_count = 0
    success_count = 0
    for i, mol_info in enumerate(molecules):
        smiles = mol_info["smiles"]
        print(f"[对接] {i+1}/{total_count}: {smiles[:40]}...")
        result = dock_molecule(smiles, center, size)
        result["smiles"] = smiles
        result.update(mol_info)
        results.append(result)
        completed_count += 1
        if result.get("success"):
            success_count += 1
        # 进度日志集成
        if hasattr(sys.stdout, "log_event"):
            sys.stdout.log_event(
                "docking_progress",
                current=completed_count,
                total=total_count,
                success_rate=success_count / max(completed_count, 1),
            )
    return results
