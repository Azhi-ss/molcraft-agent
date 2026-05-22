"""分子对接模块：使用 AutoDock Vina 计算结合自由能。

改进点（H029）:
- 引入多构象对接增强（multi-conformer docking）
- 文献依据:
  - Coscientist (Boiko et al., 2023): "performing experiments multiple times"
    消除单次实验的随机偏差
  - GNINA Benchmarking (Molecules, 2025): 构象采样质量直接影响对接精度
  - 计算化学基本原则: 构象系综对接优于单构象对接
"""
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
        ret = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if ret != 0:
        # Final fallback: random coordinates + ETKDG
        ret = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if ret != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass

    preparator = MoleculePreparation()
    setup_list = preparator.prepare(mol)
    if not setup_list:
        return None

    pdbqt_string = PDBQTWriterLegacy.write_string(setup_list[0])[0]
    if not pdbqt_string or not pdbqt_string.strip():
        return None

    if output_path is None:
        fd, output_path = tempfile.mkstemp(suffix=".pdbqt")
        os.close(fd)

    with open(output_path, "w") as f:
        f.write(pdbqt_string)

    return output_path


def smiles_to_pdbqt_multi_conformer(smiles: str, n_conformers: int = 5):
    """生成多个 3D 构象并转换为 PDBQT 文件列表（H029）。

    文献依据:
    - Coscientist (Boiko et al., 2023): 重复实验消除随机偏差
    - GNINA Benchmarking (Molecules, 2025): 构象采样影响对接精度
    - 计算化学基本原则: 构象系综对接优于单构象

    实现:
    1. RDKit ETKDGv3 生成 N 个构象
    2. MMFF94 优化每个构象
    3. Meeko 转换为 PDBQT

    Args:
        smiles: 配体 SMILES 字符串
        n_conformers: 目标构象数

    Returns:
        list[str] | None: PDBQT 文件路径列表，失败返回 None
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    # 使用 ETKDGv3 生成多个构象
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0  # 自动
    params.pruneRmsThresh = 0.5  # RMSD 剪枝阈值，保持构象多样性

    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params)
    if not conf_ids:
        # 回退到单构象
        ret = AllChem.EmbedMolecule(mol, params)
        if ret != 0:
            return None
        conf_ids = [0]

    pdbqt_paths = []
    for conf_id in conf_ids:
        try:
            # 优化当前构象
            AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, maxIters=200)
        except Exception:
            pass

        # 提取单个构象为独立 Mol 对象
        conf_mol = Chem.Mol(mol, False, conf_id)

        # 用 Meeko 准备 PDBQT
        preparator = MoleculePreparation()
        try:
            setup_list = preparator.prepare(conf_mol)
            if not setup_list:
                continue
            pdbqt_string = PDBQTWriterLegacy.write_string(setup_list[0])[0]
        except Exception:
            continue

        fd, output_path = tempfile.mkstemp(suffix=".pdbqt")
        with os.fdopen(fd, "w") as f:
            f.write(pdbqt_string)
        pdbqt_paths.append(output_path)

    return pdbqt_paths if pdbqt_paths else None


def dock_molecule(smiles: str, center=None, size=None, exhaustiveness=None, seed: int = 0,
                  n_conformers: int = 1):
    """对单个分子进行对接，返回结合能（kcal/mol）。

    H029 改进：支持多构象对接增强。
    当 n_conformers > 1 时，生成多个 3D 构象，对每个独立对接，
    取最优结合能作为结果。这模拟了"构象系综对接"策略，
    减少单构象可能错过最优结合模式的风险。

    参数:
        smiles: 配体 SMILES
        center, size: 对接盒子参数
        exhaustiveness: Vina 搜索精度
        seed: 随机种子（0 = 使用默认 42）
        n_conformers: 生成的构象数量（默认 1，即单构象；设为 >1 启用多构象）

    返回字典包含:
        - binding_energy: 最佳构象能量（负值=越好）
        - poses: 构象能量列表
        - success: 是否成功
        - error: 失败原因
        - n_conformers_tried: 实际尝试的构象数（H029 新增）
        - best_conformer_idx: 最优构象索引（H029 新增）
    """
    receptor = prepare_receptor()

    if center is None:
        center = DOCKING_CENTER
    if size is None:
        size = DOCKING_SIZE
    if exhaustiveness is None:
        exhaustiveness = DOCKING_EXHAUSTIVENESS

    vina_seed = seed if seed != 0 else 42

    # H029: 多构象对接增强
    if n_conformers > 1:
        pdbqt_paths = smiles_to_pdbqt_multi_conformer(smiles, n_conformers=n_conformers)
        if not pdbqt_paths:
            return {"success": False, "error": "多构象配体准备失败"}

        best_energy = float('inf')
        best_poses = None
        best_idx = -1
        all_successful = []

        for i, ligand_pdbqt in enumerate(pdbqt_paths):
            if ligand_pdbqt is None:
                continue
            try:
                v = Vina(sf_name="vina", seed=vina_seed, verbosity=0)
                v.set_receptor(receptor)
                v.set_ligand_from_file(ligand_pdbqt)
                v.compute_vina_maps(center=center, box_size=size)
                v.dock(exhaustiveness=exhaustiveness, n_poses=5)
                energies = v.energies(n_poses=1)
                energy = float(energies[0][0])
                all_successful.append((i, energy))
                if energy < best_energy:
                    best_energy = energy
                    best_poses = energies.tolist()
                    best_idx = i
            except Exception:
                continue
            finally:
                if ligand_pdbqt and ligand_pdbqt.startswith(tempfile.gettempdir()):
                    try:
                        os.remove(ligand_pdbqt)
                    except Exception:
                        pass

        if best_idx < 0:
            return {"success": False, "error": "所有构象对接均失败"}

        return {
            "success": True,
            "binding_energy": round(best_energy, 3),
            "poses": best_poses,
            "n_conformers_tried": len(all_successful),
            "best_conformer_idx": best_idx,
            "conformer_energies": [e for _, e in all_successful],
        }

    # 单构象模式（原有逻辑，保持向后兼容）
    ligand_pdbqt = smiles_to_pdbqt(smiles)
    if ligand_pdbqt is None:
        return {"success": False, "error": "配体准备失败"}

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
            "n_conformers_tried": 1,
            "best_conformer_idx": 0,
        }
    except Exception as e:
        if ligand_pdbqt and ligand_pdbqt.startswith(tempfile.gettempdir()):
            try:
                os.remove(ligand_pdbqt)
            except Exception:
                pass
        return {"success": False, "error": str(e)}


def dock_molecule_consensus(smiles: str, center=None, size=None, exhaustiveness=None,
                            n_runs: int = 3, seeds: list = None, n_conformers: int = 1):
    """共识对接：多次独立对接取中位数结合能（H009）。

    Coscientist 核心模式："performing experiments multiple times" —
    通过不同随机种子执行 N 次独立对接，消除单次运行的随机偏差，
    取中位数作为共识评分。

    H029 改进：支持多构象对接，每个 seed 下取多构象最优。

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
        result = dock_molecule(smiles, center, size, exhaustiveness, seed=seed,
                              n_conformers=n_conformers)
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


def batch_dock(molecules, center=None, size=None, n_conformers: int = 1):
    """批量对接分子并返回结果。

    H029 改进：支持 n_conformers 参数传递到 dock_molecule。
    """
    results = []
    total_count = len(molecules)
    completed_count = 0
    success_count = 0
    for i, mol_info in enumerate(molecules):
        smiles = mol_info["smiles"]
        # Skip multi-fragment molecules (meeko can't handle them)
        mol_check = Chem.MolFromSmiles(smiles)
        if mol_check is None or len(Chem.GetMolFrags(mol_check)) > 1:
            print(f"[对接] {i+1}/{total_count}: {smiles[:40]}... SKIP (multi-fragment)")
            result = {"smiles": smiles, "success": False, "error": "multi-fragment"}
            result.update(mol_info)
            results.append(result)
            completed_count += 1
            continue
        print(f"[对接] {i+1}/{total_count}: {smiles[:40]}...")
        try:
            result = dock_molecule(smiles, center, size, n_conformers=n_conformers)
        except Exception as e:
            print(f"  对接失败: {e}")
            result = {"smiles": smiles, "success": False, "error": str(e)}
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
