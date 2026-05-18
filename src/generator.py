"""分子生成模块：使用 RDKit 进行基于变异的分子生成。

改进点（H002）:
- 引入对接引导生成（docking-guided generation）
- 文献依据:
  - MOOSE-Chem (Yang et al., 2025): 进化算法导航组合空间
  - Coscientist (Boiko et al., 2023): 基于实验结果的迭代反思
  - 综述第4.2节: Post-Execution Feedback 策略

改进点（H010）:
- 引入 Scaffold Hopping（骨架替换）变异算子
- 使用 BRICS 分解 + BM 骨架识别 → 骨架库替换 → 侧链重连
- 文献依据:
  - Deep Lead Optimization (JACS 2024): 明确定义 Scaffold Hopping 为四个核心
    先导化合物优化子任务之一，提出"替换核心骨架同时保留有利取代基"
  - MOOSE-Chem (Yang et al., 2025): "Diverse initial population is essential
    for evolutionary search to avoid premature convergence"
  - 综述第3.2节: Scaffold hopping 是药物发现的核心策略之一
"""
import random
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Descriptors, QED, BRICS
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.DataStructs import TanimotoSimilarity
from evaluator import evaluate_molecule, passes_filters

# 压制 RDKit 错误输出
rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')


# 药物样骨架和片段库（H007: 扩充至55个，覆盖更多含氮稠环和饱和杂环）
# 文献依据:
#   - MOOSE-Chem (Yang et al., 2025): "Diverse initial population is essential
#     for evolutionary search to avoid premature convergence"
#   - ChemCrow (Bran et al., 2024): 工具/库的丰富度直接决定 Agent 探索的化学空间边界
#   - 综述第3.2节: Scaffold hopping 是药物发现的核心策略之一
SCAFFOLDS = [
    # ===== 单环芳烃（6个）=====
    "c1ccc(cc1)",           # 苯环
    "c1ccccc1C",            # 甲苯
    "c1ccc(cc1)O",          # 苯酚
    "c1ccc(cc1)N",          # 苯胺
    "COc1ccc(cc1)",         # 苯甲醚
    "c1ccc(cc1)CN",         # 苄胺

    # ===== 多环芳烃（2个）=====
    "c1ccc2ccccc2c1",       # 萘
    "c1ccc2c(c1)cccc2",     # 萘变体

    # ===== 含氮单杂环（4个）=====
    "c1ccc(nc1)",           # 吡啶
    "c1cncnc1",             # 嘧啶
    "c1cnccc1",             # 吡啶变体
    "c1c[nH]cn1",           # 咪唑

    # ===== 饱和杂环（8个）=====
    "C1CCOC1",              # 四氢呋喃
    "C1CCNC1",              # 吡咯烷
    "C1CCNCC1",             # 哌啶
    "C1CNCCN1",             # 哌嗪
    "C1COCCN1",             # 吗啉
    "C1CSCN1",              # 硫代吗啉
    "C1COC1",               # 氧杂环丁烷
    "C1CNC1",               # 氮杂环丁烷

    # ===== 含氮稠环（12个）=====
    "c1ccc2c(c1)Ncnc2",     # 喹唑啉
    "c1ccc2c(c1)ncn2",      # 苯并咪唑
    "c1ccc2c(c1)cccn2",     # 吲哚
    "c1ccc2c(c1)ocn2",      # 苯并噁唑
    "c1ccc2c(c1)scn2",      # 苯并噻唑
    "c1ccc2ncccc2c1",       # 喹啉
    "c1ccc2ccncc2c1",       # 异喹啉
    "c1ccc2c(c1)[nH]c2",    # 吲哚啉（二氢吲哚）
    "c1ccc2c(c1)CCN2",      # 二氢吲哚（含N）
    "c1ccc2c(c1)CCNC2",     # 四氢异喹啉
    "c1ccc2c(c1)N=CN2",     # 苯并咪唑啉
    "c1nc2c([nH]1)cccc2",   # 苯并咪唑变体

    # ===== 药物常见骨架（6个）=====
    "c1ncnc2c1ncn2",        # 嘌呤
    "c1cnc2ncncc2n1",       # 蝶啶
    "c1ccc2c(c1)cncn2",     # 喹唑啉变体
    "c1ccc2c(c1)ncnc2",     # 喹唑啉（另一表示）
    "c1ccc2c(c1)OCCO2",     # 苯并二噁烷
    "c1ccc2c(c1)CCO2",      # 苯并呋喃烷

    # ===== 酰胺/磺酰胺类（6个）=====
    "c1ccc(cc1)C(=O)O",     # 苯甲酸
    "c1ccc(cc1)C(=O)N",     # 苯甲酰胺
    "c1ccc(cc1)S(=O)(=O)N", # 磺酰胺
    "c1ccc(cc1)C(=O)Nc2ccccc2",  # 二苯甲酮酰胺
    "c1ccc(cc1)NC(=O)c2ccccc2",  # N-苯基苯甲酰胺
    "Cc1ccc(cc1)S(=O)(=O)Nc2ccccc2",  # 对甲苯磺酰苯胺

    # ===== 卤代/其他（5个）=====
    "c1cc(ccc1F)F",         # 二氟苯
    "c1cc(ccc1Cl)Cl",       # 二氯苯
    "c1ccc(cc1)CC(=O)O",    # 苯乙酸
    "c1ccc(cc1)C(=O)c2ccccc2",  # 二苯甲酮
    "c1ccc(cc1)OCc2ccccc2",     # 二苯甲醚

    # ===== 饱和环-芳环稠合（6个）=====
    "c1ccc2c(c1)CCCC2",     # 四氢萘
    "c1ccc2c(c1)CCCCC2",    # 十氢萘骨架
    "c1ccc2c(c1)NCCC2",     # 四氢喹啉
    "c1ccc2c(c1)OCC2",      # 2,3-二氢苯并呋喃
    "c1ccc2c(c1)SCC2",      # 2,3-二氢苯并噻吩
    "c1ccc2c(c1)CCO2",      # 苯并二氢吡喃
]

LINKERS = [
    "",                     # 直接连接
    "C",                    # 亚甲基
    "CC",                   # 亚乙基
    "O",                    # 醚键
    "NH",                   # 胺键
    "C(=O)",               # 羰基
    "C(=O)N",              # 酰胺
    "C(=O)O",              # 酯键
    "S(=O)(=O)",           # 砜
    "S(=O)(=O)N",          # 磺酰胺
    "NHC(=O)",             # 反向酰胺
    "OC(=O)",              # 反向酯键
    "C#C",                 # 炔键
    "C=C",                 # 烯键
]


def random_mutate_smiles(smiles: str, n_mutations: int = 1):
    """对 SMILES 字符串应用随机变异。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    for _ in range(n_mutations):
        mol = _mutate_mol(mol)
        if mol is None:
            return None

    new_smiles = Chem.MolToSmiles(mol, canonical=True)
    return new_smiles


def _mutate_mol(mol):
    """单次变异：添加/替换/删除/连接/骨架替换。

    H010 改进：新增 scaffold hopping 算子（~20% 概率），
    对应 Deep Lead Optimization 的四个核心子任务中的 Scaffold Hopping。
    """
    choice = random.random()
    try:
        if choice < 0.25:
            mol = _add_substituent(mol)
        elif choice < 0.50:
            mol = _replace_atom(mol)
        elif choice < 0.60:
            mol = _remove_terminal(mol)
        elif choice < 0.80:
            mol = _insert_linker(mol)
        else:
            # H010: Scaffold Hopping（20% 概率）
            result = _scaffold_hop(mol)
            if result is not None:
                mol = result
            else:
                # 骨架替换失败 → 回退到插入连接子
                mol = _insert_linker(mol)
    except Exception:
        return None
    if mol is None:
        return None
    # 返回前验证分子有效性
    try:
        Chem.SanitizeMol(mol)
        # 测试 SMILES 往返
        s = Chem.MolToSmiles(mol, canonical=True)
        m2 = Chem.MolFromSmiles(s)
        if m2 is None:
            return None
    except Exception:
        return None
    return mol


def _add_substituent(mol):
    """在随机碳原子上添加小取代基（F, Cl, OH, NH2, CH3）。"""
    substituents = ["F", "Cl", "[OH]", "[NH2]", "C"]
    subst = random.choice(substituents)
    subst_mol = Chem.MolFromSmiles(subst)
    if subst_mol is None:
        return mol

    # 查找候选原子（非氢、非末端碳原子）
    atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetDegree() < 4]
    if not atoms:
        return mol
    atom = random.choice(atoms)

    # 将取代基连接到分子上
    combo = Chem.CombineMols(mol, subst_mol)
    emol = Chem.EditableMol(combo)
    new_bond_idx = combo.GetNumAtoms() - 1
    emol.AddBond(atom.GetIdx(), new_bond_idx, Chem.BondType.SINGLE)
    new_mol = emol.GetMol()
    Chem.SanitizeMol(new_mol)
    return new_mol


def _replace_atom(mol):
    """随机替换一个杂原子或碳原子为另一种原子。"""
    replacements = {6: [7, 8, 9, 16], 7: [6, 8], 8: [6, 7, 16], 16: [6, 8]}
    atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() in replacements]
    if not atoms:
        return mol
    atom = random.choice(atoms)
    new_atom_num = random.choice(replacements[atom.GetAtomicNum()])
    atom.SetAtomicNum(new_atom_num)
    Chem.SanitizeMol(mol)
    return mol


def _remove_terminal(mol):
    """删除一个末端原子。"""
    terminals = [a for a in mol.GetAtoms() if a.GetDegree() == 1 and not a.IsInRing()]
    if not terminals:
        return mol
    atom = random.choice(terminals)
    emol = Chem.EditableMol(mol)
    emol.RemoveAtom(atom.GetIdx())
    new_mol = emol.GetMol()
    Chem.SanitizeMol(new_mol)
    return new_mol


def _insert_linker(mol):
    """在两个片段之间插入一个小连接子。"""
    linkers = ["C", "O", "N"]
    linker = random.choice(linkers)
    linker_mol = Chem.MolFromSmiles(linker)
    if linker_mol is None:
        return mol

    combo = Chem.CombineMols(mol, linker_mol)
    emol = Chem.EditableMol(combo)
    atoms = list(mol.GetAtoms())
    if len(atoms) < 2:
        return mol
    a1, a2 = random.sample(atoms, 2)
    linker_idx = combo.GetNumAtoms() - 1
    emol.AddBond(a1.GetIdx(), linker_idx, Chem.BondType.SINGLE)
    emol.AddBond(a2.GetIdx(), linker_idx, Chem.BondType.SINGLE)
    new_mol = emol.GetMol()
    try:
        Chem.SanitizeMol(new_mol)
    except Exception:
        return mol
    return new_mol


def _scaffold_hop(mol):
    """Scaffold Hopping 算子（H010）：替换核心骨架，保留侧链装饰。

    文献依据：
    - Deep Lead Optimization (JACS 2024, §3.1): Scaffold Hopping 是先导化合物
      优化的四个核心子任务之一，定义为"替换核心骨架同时保留有利取代基"
    - 综述 §3.2: Scaffold hopping 是命中化合物到先导化合物过渡中的关键操作

    算法流程：
    1. BRICS 分解分子 → 得到片段列表
    2. 识别最大片段作为"核心骨架"
    3. 从 SCAFFOLDS 库选择不同于当前骨架的新骨架
    4. 在新骨架上寻找合理连接点，将侧链重连
    5. 验证生成分子的化学合理性
    """
    if mol is None or mol.GetNumAtoms() < 8:
        # 过小的分子无法做有意义的骨架替换
        return None

    try:
        # Step 1: BRICS 分解（使用 BRICSDecompose 返回 SMILES 集合）
        from rdkit.Chem.BRICS import BRICSDecompose
        frag_smiles_set = BRICSDecompose(mol)

        if not frag_smiles_set or len(frag_smiles_set) < 2:
            # BRICS 无法分解 → 回退到 Murcko 骨架分析
            core = MurckoScaffold.GetScaffoldForMol(mol)
            if core is None or core.GetNumAtoms() < 3:
                return None
            # 用 Murcko 骨架作为核心
            frag_list = [core]
            # 获取侧链：移除 Murcko 骨架后剩余部分
            side_chains = Chem.DeleteSubstructs(mol, core)
            if side_chains is not None and side_chains.GetNumAtoms() > 0:
                frag_list.append(side_chains)
            if len(frag_list) < 2:
                return None
            # Murcko 路径：直接用 mol 对象，跳过 BRICS 处理
            valid_frags = frag_list
            side_frags_raw = [frag_list[1]] if len(frag_list) > 1 else []
            core_frag = frag_list[0]
        else:
            # BRICS 路径：将 SMILES 转换为 Mol 对象
            frag_mols = {}
            for smi in frag_smiles_set:
                fm = Chem.MolFromSmiles(smi)
                if fm is not None and fm.GetNumAtoms() >= 2:
                    # 计算非 dummy 原子数
                    nd = sum(1 for a in fm.GetAtoms() if a.GetAtomicNum() != 0)
                    if nd >= 2:
                        frag_mols[fm] = nd

            if len(frag_mols) < 2:
                return None

            # Step 2: 按非 dummy 原子数排序，最大的作为核心骨架
            sorted_frags = sorted(frag_mols.items(), key=lambda x: x[1], reverse=True)
            core_frag = sorted_frags[0][0]

            # 侧链：保留原始 mol（含 dummy atom）用于识别连接点
            side_frags_raw = [f for f, _ in sorted_frags[1:]]

            empty = set()
            valid_frags = [f for f, _ in sorted_frags]
            # Make core_smiles for later comparison
            core_smiles_nodummy = Chem.MolToSmiles(core_frag, canonical=True)

        # Get core SMILES (for similarity check later)
        try:
            core_smiles = Chem.MolToSmiles(core_frag, canonical=True)
        except Exception:
            core_smiles = ""

        # Step 3: 从 SCAFFOLDS 库中选择新骨架
        # 排除与当前骨架过于相似的（Tanimoto < 0.4）以鼓励多样性
        candidate_scaffolds = []
        core_fp = None
        try:
            core_fp = AllChem.GetMorganFingerprintAsBitVect(core_frag, 2, nBits=1024)
        except Exception:
            pass

        for s_smi in SCAFFOLDS:
            if s_smi == core_smiles:
                continue
            if core_fp is not None:
                try:
                    s_mol = Chem.MolFromSmiles(s_smi)
                    if s_mol is None:
                        continue
                    s_fp = AllChem.GetMorganFingerprintAsBitVect(s_mol, 2, nBits=1024)
                    sim = TanimotoSimilarity(core_fp, s_fp)
                    if sim < 0.4:  # 选择足够不同的骨架
                        candidate_scaffolds.append(s_smi)
                except Exception:
                    candidate_scaffolds.append(s_smi)
            else:
                candidate_scaffolds.append(s_smi)

        if not candidate_scaffolds:
            # 如果没有足够不同的骨架，允许更相似的
            candidate_scaffolds = [s for s in SCAFFOLDS if s != core_smiles]

        if not candidate_scaffolds:
            return None

        new_scaffold_smi = random.choice(candidate_scaffolds)
        new_scaffold = Chem.MolFromSmiles(new_scaffold_smi)
        if new_scaffold is None:
            return None

        # Step 4: 将侧链连接到新骨架
        # 在新骨架上找可连接原子（碳原子，度 < 4）
        scaffold_atoms = []
        for atom in new_scaffold.GetAtoms():
            if atom.GetAtomicNum() == 6 and atom.GetDegree() < 4:
                scaffold_atoms.append(atom)

        if not scaffold_atoms:
            # 放宽到任何度 < 4 的原子
            for atom in new_scaffold.GetAtoms():
                if atom.GetDegree() < 4 and atom.GetAtomicNum() in (6, 7, 8):
                    scaffold_atoms.append(atom)

        if not scaffold_atoms:
            return None

        # 逐个连接侧链（最多连接3个，避免过度复杂）
        combo = new_scaffold
        n_sides = min(len(side_frags_raw), 3)
        attached_atoms = set()

        for i in range(n_sides):
            side = side_frags_raw[i]

            # Step 4a: 找到侧链的连接点（与 dummy atom 相邻的真实原子）
            side_attachment_atoms = []
            dummy_indices = []
            for atom in side.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    dummy_indices.append(atom.GetIdx())
                    for nbr in atom.GetNeighbors():
                        if nbr.GetAtomicNum() != 0:
                            side_attachment_atoms.append(nbr)

            if not side_attachment_atoms:
                # 没有 dummy atom 时，退回到通用连接点查找
                for atom in side.GetAtoms():
                    if atom.GetAtomicNum() != 0 and atom.GetDegree() < 4:
                        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
                            side_attachment_atoms.insert(0, atom)
                        elif atom.GetAtomicNum() in (6, 7, 8) and atom.GetDegree() < 4:
                            side_attachment_atoms.append(atom)

            if not side_attachment_atoms:
                continue

            # Step 4b: 从侧链中移除 dummy atom，得到干净的侧链片段
            side_clean = Chem.RWMol(side)
            for di in sorted(dummy_indices, reverse=True):
                side_clean.RemoveAtom(di)
            side_clean = side_clean.GetMol()
            # 不进行 sanitize，保留自由基用于成键

            # 由于移除了 dummy 原子，侧链中原子索引可能偏移
            # 重新计算连接点在干净侧链中的索引
            old_to_new = {}
            new_idx = 0
            for old_idx in range(side.GetNumAtoms()):
                if old_idx not in dummy_indices:
                    old_to_new[old_idx] = new_idx
                    new_idx += 1

            side_attach_new_indices = [
                old_to_new[a.GetIdx()]
                for a in side_attachment_atoms
                if a.GetIdx() in old_to_new
            ]
            if not side_attach_new_indices:
                continue

            # Step 4c: 选择骨架上的可用连接点
            available_scaffold = [a for a in scaffold_atoms
                                  if a.GetIdx() not in attached_atoms]
            if not available_scaffold:
                break

            scaffold_anchor = random.choice(available_scaffold)
            side_anchor_new_idx = random.choice(side_attach_new_indices)

            # Step 4d: 组合分子并添加键
            new_combo = Chem.CombineMols(combo, side_clean)
            emol = Chem.EditableMol(new_combo)

            combo_natoms = combo.GetNumAtoms()
            side_anchor_global_idx = combo_natoms + side_anchor_new_idx

            emol.AddBond(scaffold_anchor.GetIdx(), side_anchor_global_idx,
                         Chem.BondType.SINGLE)

            try:
                combo = emol.GetMol()
                Chem.SanitizeMol(combo)
                attached_atoms.add(scaffold_anchor.GetIdx())
            except Exception:
                # 这个侧链连接失败，继续尝试下一个
                continue

        # Step 5: 最终验证
        try:
            Chem.SanitizeMol(combo)
            new_smiles = Chem.MolToSmiles(combo, canonical=True)
            # SMILES 往返验证
            verify_mol = Chem.MolFromSmiles(new_smiles)
            if verify_mol is None:
                return None
            # 确保生成的是与原始分子不同的结构
            orig_smiles = Chem.MolToSmiles(mol, canonical=True)
            if new_smiles == orig_smiles:
                return None
            return combo
        except Exception:
            return None

    except Exception:
        return None

def generate_molecules(strategy="mutate", n_molecules=50, scaffold=None):
    """生成候选药物分子。

    策略:
        - "mutate": 从种子骨架变异
        - "combine": 用连接子组合骨架
        - "random": 随机 SMILES 生成（非常基础）
    """
    molecules = set()
    attempts = 0
    max_attempts = n_molecules * 20

    if strategy == "mutate":
        seeds = SCAFFOLDS if scaffold is None else [scaffold]
        while len(molecules) < n_molecules and attempts < max_attempts:
            attempts += 1
            seed = random.choice(seeds)
            n_mut = random.randint(1, 4)
            new_smiles = random_mutate_smiles(seed, n_mut)
            if new_smiles and new_smiles not in molecules:
                props = evaluate_molecule(new_smiles)
                if passes_filters(props):
                    molecules.add(new_smiles)

    elif strategy == "combine":
        while len(molecules) < n_molecules and attempts < max_attempts:
            attempts += 1
            n_frag = random.randint(2, 3)
            frags = random.sample(SCAFFOLDS, n_frag)
            linkers = random.sample(LINKERS, n_frag - 1)

            # 用连接子将片段拼接成 SMILES
            parts = []
            for i, frag in enumerate(frags):
                parts.append(frag)
                if i < len(linkers):
                    parts.append(linkers[i])
            combined = "".join(parts)

            mol = Chem.MolFromSmiles(combined)
            if mol is not None:
                smiles = Chem.MolToSmiles(mol, canonical=True)
                if smiles not in molecules:
                    props = evaluate_molecule(smiles)
                    if passes_filters(props):
                        molecules.add(smiles)

    elif strategy == "random":
        # 非常基础：随机组合片段
        while len(molecules) < n_molecules and attempts < max_attempts:
            attempts += 1
            frag = random.choice(SCAFFOLDS)
            new_smiles = random_mutate_smiles(frag, random.randint(2, 5))
            if new_smiles and new_smiles not in molecules:
                props = evaluate_molecule(new_smiles)
                if passes_filters(props):
                    molecules.add(new_smiles)

    result = [evaluate_molecule(s) for s in molecules]
    result.sort(key=lambda x: x["qed"], reverse=True)
    return result


def _diverse_selection(candidates, n_select, diversity_weight=0.3):
    """多样性保持选择（H011）：使用贪心 MMD 算法选择候选分子。

    文献依据:
    - MOOSE-Chem (Yang et al., 2025): "Diverse initial population is essential
      for evolutionary search to avoid premature convergence"
    - MolLEO (Wang et al., 2024b): 多目标进化优化中多样性是核心维度

    算法: Maximum Minimal Distance (MMD) 贪心选择
    1. 计算所有分子的 Morgan 指纹
    2. 从结合能最优分子开始
    3. 迭代选择: 综合评分 = (1-w)×BE_norm + w×diversity_norm
       其中 diversity_norm = min(Tanimoto_distance to already selected)

    Args:
        candidates: list[dict], 每个含 "smiles", "binding_energy"
        n_select: 选择数量
        diversity_weight: 多样性权重 (0=纯BE, 1=纯多样性)

    Returns:
        list[dict]: 按综合评分排序的选择结果
    """
    if len(candidates) <= n_select:
        return list(candidates)

    # 预计算指纹
    fps = []
    for c in candidates:
        try:
            mol = Chem.MolFromSmiles(c["smiles"])
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            fps.append(fp)
        except Exception:
            fps.append(None)

    # 归一化结合能到 [0, 1]（越低越好 → 归一化后 1 为最好）
    energies = [c.get("binding_energy", 999) for c in candidates]
    e_min, e_max = min(energies), max(energies)
    if e_max > e_min:
        be_norm = [(e_max - e) / (e_max - e_min) for e in energies]
    else:
        be_norm = [0.5] * len(energies)

    selected = []
    selected_indices = set()

    # Step 1: 选择结合能最优分子作为起点
    best_idx = min(range(len(candidates)), key=lambda i: energies[i])
    selected.append(candidates[best_idx])
    selected_indices.add(best_idx)

    # Step 2: 贪心选择剩余分子
    for _ in range(1, n_select):
        best_score = -1
        best_idx = -1

        for i, cand in enumerate(candidates):
            if i in selected_indices:
                continue
            if fps[i] is None:
                continue

            # 计算与已选集合的最小 Tanimoto 距离（距离 = 1 - 相似度）
            max_sim = 0.0
            for si in selected_indices:
                if fps[si] is not None:
                    sim = TanimotoSimilarity(fps[i], fps[si])
                    if sim > max_sim:
                        max_sim = sim
            diversity_norm = 1.0 - max_sim  # 距离作为多样性分数

            # 综合评分
            score = (1 - diversity_weight) * be_norm[i] + diversity_weight * diversity_norm
            if score > best_score:
                best_score = score
                best_idx = i

        if best_idx >= 0:
            selected.append(candidates[best_idx])
            selected_indices.add(best_idx)

    # 按结合能排序返回
    selected.sort(key=lambda x: x.get("binding_energy", 999))
    return selected


def generate_with_docking_guidance(
    docking_fn,
    n_molecules=50,
    batch_size=10,
    n_generations=3,
    top_k=5,
    strategy="mutate",
    scaffold=None,
):
    """对接引导的分子生成（H002 + H011 多样性保持）。

    核心思想：将分子对接作为适应度函数，嵌入生成循环中。
    每批生成少量分子 → 对接评估 → 多样性保持选择 → 变异产生下一代。
    这避免了盲生成大量低质量分子，同时防止过早收敛。

    H011 改进：种子选择加入多样性保持（贪心 MMD 算法），
    避免纯结合能选择导致的过早收敛。

    Args:
        docking_fn: 对接函数，接收 SMILES 字符串，返回 dict 包含 "binding_energy"
        n_molecules: 目标生成分子总数
        batch_size: 每批生成的候选分子数
        n_generations: 进化代数
        top_k: 每代按结合能保留的 top_k，剩余 (n_seeds - top_k) 由多样性选择
        strategy: 初始生成策略
        scaffold: 可选种子骨架

    Returns:
        list[dict]: 通过过滤且对接成功的分子信息列表，按结合能排序
    """
    import sys

    all_evaluated = []  # 所有经过对接评估的分子
    current_seeds = None

    for gen in range(n_generations):
        batch_mols = []
        attempts = 0
        max_attempts = batch_size * 30

        if gen == 0:
            # 初始代：用传统策略生成
            seeds = SCAFFOLDS if scaffold is None else [scaffold]
            while len(batch_mols) < batch_size and attempts < max_attempts:
                attempts += 1
                seed = random.choice(seeds)
                n_mut = random.randint(1, 4)
                new_smiles = random_mutate_smiles(seed, n_mut)
                if new_smiles and new_smiles not in {m["smiles"] for m in batch_mols}:
                    props = evaluate_molecule(new_smiles)
                    if passes_filters(props):
                        batch_mols.append(props)
        else:
            # 后续代：从种子变异
            seed_smiles_list = [s["smiles"] for s in current_seeds]
            while len(batch_mols) < batch_size and attempts < max_attempts:
                attempts += 1
                seed_smiles = random.choice(seed_smiles_list)
                n_mut = random.randint(1, 3)  # 后续代变异强度略低，保持稳定性
                new_smiles = random_mutate_smiles(seed_smiles, n_mut)
                if new_smiles and new_smiles not in {m["smiles"] for m in batch_mols}:
                    props = evaluate_molecule(new_smiles)
                    if passes_filters(props):
                        batch_mols.append(props)

        # 对接评估
        docked_batch = []
        for mol_info in batch_mols:
            smiles = mol_info["smiles"]
            try:
                dock_result = docking_fn(smiles)
                if dock_result.get("success"):
                    mol_info["binding_energy"] = dock_result.get("binding_energy")
                    mol_info["docking_success"] = True
                    docked_batch.append(mol_info)
            except Exception:
                continue

        # 按结合能排序（越低越好）
        docked_batch.sort(key=lambda x: x.get("binding_energy", 999))

        # 累积到全局池
        all_evaluated.extend(docked_batch)

        # H011: 多样性保持选择种子
        # 策略: 50% 种子来自纯 BE（确保选择压力），50% 来自多样性选择（确保探索性）
        n_be_seeds = max(1, top_k // 2)
        n_diverse_seeds = max(1, top_k - n_be_seeds)
        n_total_seeds = min(top_k, len(docked_batch))

        # 纯 BE 种子
        be_seeds = docked_batch[:n_be_seeds]

        # 多样性种子（从剩余分子中选择）
        remaining = docked_batch[n_be_seeds:]
        if remaining and n_diverse_seeds > 0:
            diverse_seeds = _diverse_selection(
                remaining, min(n_diverse_seeds, len(remaining)), diversity_weight=0.5
            )
        else:
            diverse_seeds = []

        # 合并种子（去重）
        seed_smiles_set = {s["smiles"] for s in be_seeds}
        current_seeds = list(be_seeds)
        for s in diverse_seeds:
            if s["smiles"] not in seed_smiles_set:
                seed_smiles_set.add(s["smiles"])
                current_seeds.append(s)

        # 确保有足够种子
        if len(current_seeds) < n_total_seeds:
            for cand in docked_batch:
                if cand["smiles"] not in seed_smiles_set:
                    seed_smiles_set.add(cand["smiles"])
                    current_seeds.append(cand)
                if len(current_seeds) >= n_total_seeds:
                    break

        current_seeds = current_seeds[:n_total_seeds]

        # 计算多样性指标
        if len(docked_batch) >= 2:
            avg_sim = _avg_pairwise_similarity(docked_batch[:min(20, len(docked_batch))])
        else:
            avg_sim = 0.0

        print(
            f"[H002+H011] Gen {gen+1}/{n_generations}: generated {len(batch_mols)}, "
            f"docked {len(docked_batch)}, best {docked_batch[0]['binding_energy'] if docked_batch else 'N/A'}, "
            f"seeds {len(current_seeds)} (BE={n_be_seeds} diverse={len(diverse_seeds)}), "
            f"avg_pairwise_sim={avg_sim:.3f}",
            file=sys.stderr,
        )

    # 全局去重 + 多样性最终排序
    seen = set()
    unique = []
    for m in all_evaluated:
        smi = m["smiles"]
        if smi not in seen:
            seen.add(smi)
            unique.append(m)

    unique.sort(key=lambda x: x.get("binding_energy", 999))

    # H011: 最终选择也加入多样性 — 从 top 2*n_molecules 中做多样性选择
    # 确保最终的分子集合既有高结合能又有多样性
    candidate_pool = unique[:min(n_molecules * 2, len(unique))]
    if len(candidate_pool) > n_molecules:
        # 保留 top 60% 纯BE + 40% 多样性
        n_be = int(n_molecules * 0.6)
        n_div = n_molecules - n_be
        be_final = candidate_pool[:n_be]
        remaining_pool = candidate_pool[n_be:]
        if remaining_pool and n_div > 0:
            div_final = _diverse_selection(
                remaining_pool, min(n_div, len(remaining_pool)), diversity_weight=0.5
            )
        else:
            div_final = []

        # 合并: BE 种子 + 多样性种子
        final_set = set(m["smiles"] for m in be_final)
        result = list(be_final)
        for m in div_final:
            if m["smiles"] not in final_set:
                final_set.add(m["smiles"])
                result.append(m)
        # 如果不够，从剩余中补充
        for m in candidate_pool:
            if len(result) >= n_molecules:
                break
            if m["smiles"] not in final_set:
                final_set.add(m["smiles"])
                result.append(m)
        result.sort(key=lambda x: x.get("binding_energy", 999))
        return result[:n_molecules]

    return candidate_pool[:n_molecules]


def _avg_pairwise_similarity(molecules):
    """计算分子集合的平均成对 Tanimoto 相似度（H011 辅助）。"""
    if len(molecules) < 2:
        return 0.0
    fps = []
    for m in molecules:
        try:
            mol = Chem.MolFromSmiles(m.get("smiles", ""))
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            fps.append(fp)
        except Exception:
            fps.append(None)

    sims = []
    for i in range(len(fps)):
        for j in range(i + 1, len(fps)):
            if fps[i] is not None and fps[j] is not None:
                sims.append(TanimotoSimilarity(fps[i], fps[j]))
    return sum(sims) / len(sims) if sims else 0.0
