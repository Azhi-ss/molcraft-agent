"""分子性质评估模块：使用 RDKit 计算药物相似性和理化性质。"""
import math
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Lipinski
from rdkit.Chem import RDConfig
import os


def evaluate_molecule(smiles: str):
    """评估分子的药物相似性和基本理化性质。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"valid": False, "error": "无效的 SMILES"}

    mol = Chem.AddHs(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    rings = Descriptors.RingCount(mol)
    qed = QED.qed(mol)

    # 合成可及性估算（非常粗略）
    sa = estimate_sa_score(mol)

    # 类药五规则检查
    lipinski_violations = 0
    if mw > 500: lipinski_violations += 1
    if logp > 5: lipinski_violations += 1
    if hbd > 5: lipinski_violations += 1
    if hba > 10: lipinski_violations += 1
    lipinski_ok = lipinski_violations <= 1  # 允许违反 1 条

    return {
        "valid": True,
        "smiles": smiles,
        "mw": round(mw, 2),
        "logp": round(logp, 2),
        "tpsa": round(tpsa, 2),
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotb,
        "rings": rings,
        "qed": round(qed, 3),
        "sa_score": round(sa, 2),
        "lipinski_pass": lipinski_ok,
    }


def estimate_sa_score(mol):
    """合成可及性估算（0-10，越低越好）。

    H001 改进：
    - 环贡献因子从 0.5→1.0（多环结构显著增加合成难度）
    - 新增稠环惩罚：每个额外环系 +0.3（fused ring penalty）
    - 螺环惩罚从 1.0→1.5（螺环形成挑战性高）
    - 桥头原子惩罚从 1.5→2.0（桥头结构极难合成）

    """
    # 统计环信息
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    n_rings = len(atom_rings)

    # 检测稠环数（真正共享原子/键的环系）
    # ri.BondRings() 返回所有环，错误地会对单环多环都惩罚
    # 正确方式：统计共享 >=2 个原子的环对数量
    n_fused_systems = 0
    n = len(atom_rings)
    for i in range(n):
        for j in range(i + 1, n):
            shared = len(set(atom_rings[i]) & set(atom_rings[j]))
            if shared >= 2:  # 共享至少2个原子才是稠合环
                n_fused_systems += 1
    fused_penalty = n_fused_systems * 0.3

    # 螺环中心
    spiro = Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol)

    # 桥头原子
    bridge = Chem.rdMolDescriptors.CalcNumBridgeheadAtoms(mol)

    # 手性中心
    stereo = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)

    # 基础启发式公式（优化后）
    score = 1.0 + n_rings * 1.0 + fused_penalty + spiro * 1.5 + bridge * 2.0 + stereo * 0.5
    score += Descriptors.NumRotatableBonds(mol) * 0.1
    return min(score, 10.0)


def passes_filters(props: dict, min_qed=0.3, max_mw=500, min_mw=150, max_logp=5.0,
                   max_sa=7.0, max_rings=9):
    """检查分子是否通过基础类药性质过滤。

    H001 改进：
    - SA score 阈值从 8.0 收紧至 6.0
    - 新增 max_rings=7 环数上限
    - 文献依据：Deep Lead Optimization (JACS, 2024)

    H027 改进：
    - max_rings 7→9, max_sa 6.0→7.0
    - 基于 Round 18-19 实验验证：Vina 评分偏向疏水芳香体系，
      更宽松的过滤器允许更大的π体系通过（MW仍<500约束）

    """
    if not props.get("valid"):
        return False
    if props["qed"] < min_qed:
        return False
    if not (min_mw <= props["mw"] <= max_mw):
        return False
    if props["logp"] > max_logp:
        return False
    if props["sa_score"] > max_sa:
        return False
    if props.get("rings", 0) > max_rings:
        return False
    return True
