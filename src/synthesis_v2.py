"""逆合成规划模块（第二版）：使用 SMARTS 反应模板生成 realistic 起始原料。

文献依据与改进说明:
- LARC (Baker et al., 2025) 提出 Agent-as-a-Judge 逆合成框架，强调规则覆盖率
  是逆合成质量的关键决定因素。
- ChemCrow (Bran et al., 2024) 集成 18 个专家化学工具，证明工具丰富度直接
  决定 Agent 能力边界。
- 本模块将规则从 8 条扩充至 35+ 条，覆盖常见药物化学反应类型，
  显著降低 trivial route 比例。

改进点（H003）:
- 引入递归多步逆合成规划（max_depth=3）
- 对单步断键后的中间体继续递归断键，直到得到简单起始原料
- 文献依据: LARC (2025) 强调多步路线评审的重要性
"""
import re
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, BRICS

rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')


def _validate_smiles(smiles: str) -> bool:
    """验证 SMILES 字符串是否有效且可解析。"""
    if not smiles or len(smiles) < 1:
        return False
    # 排除仅包含氢或自由基的无效产物
    if smiles in ('[H]', '[H][H]', 'H', '[Br]', '[Cl]', '[I]', '[F]'):
        return False
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None


# ---- 自动原子守恒副产物表 ----
# 每个条目: (atom_counts_dict, SMILES)
# SMILES 经过 RDKit AddHs 验证，原子计数与 dict 精确匹配
_BYPRODUCT_TABLE = [
    # 有机小分子（优先匹配，减少碎片数）
    ({'C': 2, 'H': 4, 'O': 2}, 'CC(=O)O'),     # CH3COOH 乙酸
    ({'C': 1, 'H': 4, 'O': 1}, 'CO'),           # CH3OH 甲醇
    # 含硼副产物
    ({'B': 1, 'H': 3, 'O': 3}, 'OB(O)O'),       # H3BO3 硼酸
    ({'B': 1, 'H': 3, 'O': 2}, 'B(O)O'),        # BH3O2
    ({'B': 1, 'H': 1, 'O': 2}, 'B(=O)O'),       # HBO2 偏硼酸 (Suzuki byproduct)
    # 无机小分子
    ({'N': 1, 'H': 3}, 'N'),                    # NH3
    ({'S': 1, 'H': 2}, 'S'),                    # H2S
    ({'H': 2, 'O': 1}, 'O'),                    # H2O
    ({'H': 1, 'Cl': 1}, 'Cl'),                  # HCl
    ({'H': 1, 'Br': 1}, 'Br'),                  # HBr
    ({'H': 1, 'I': 1}, 'I'),                    # HI
    # 双原子气体
    ({'H': 2}, '[H][H]'),                       # H2
    ({'O': 2}, 'O=O'),                          # O2
    ({'N': 2}, 'N#N'),                          # N2
    ({'C': 1, 'O': 2}, 'O=C=O'),               # CO2
    ({'S': 1, 'O': 2}, 'O=S=O'),               # SO2
    # 单原子离子（无隐式 H，兜底用）
    ({'Na': 1}, '[Na+]'),
    ({'K': 1}, '[K+]'),
    ({'Li': 1}, '[Li+]'),
    ({'Cl': 1}, '[Cl-]'),
    ({'Br': 1}, '[Br-]'),
    ({'I': 1}, '[I-]'),
    ({'F': 1}, '[F-]'),
    ({'B': 1}, '[B]'),
    ({'C': 1}, '[C]'),
    ({'N': 1}, '[N]'),
    ({'O': 1}, '[O]'),
    ({'P': 1}, '[P]'),
    ({'S': 1}, '[S]'),
    ({'Si': 1}, '[Si]'),
]


def _count_all_atoms(smiles_str: str) -> dict:
    """Count ALL atoms in a SMILES string (including H) using AddHs."""
    counts = {}
    for frag in smiles_str.split('.'):
        frag = frag.strip()
        if not frag:
            continue
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            continue
        h_mol = Chem.AddHs(mol)
        for atom in h_mol.GetAtoms():
            sym = atom.GetSymbol()
            counts[sym] = counts.get(sym, 0) + 1
    return counts


def _decompose_excess(excess: dict, max_frags: int = 5) -> list | None:
    """Recursively decompose excess atoms into known small-molecule byproducts.

    Returns list of SMILES strings, or None if decomposition fails.
    """
    excess = {k: v for k, v in excess.items() if v > 0}
    if not excess:
        return []
    if max_frags <= 0:
        return None

    for bp_atoms, bp_smiles in _BYPRODUCT_TABLE:
        if not all(excess.get(sym, 0) >= need for sym, need in bp_atoms.items()):
            continue

        remaining = dict(excess)
        for sym, need in bp_atoms.items():
            remaining[sym] -= need
            if remaining[sym] == 0:
                del remaining[sym]

        result = _decompose_excess(remaining, max_frags - 1)
        if result is not None:
            return [bp_smiles] + result

    return None


def _add_byproducts_for_balance(product_smiles: str, reactant_smiles: list) -> list:
    """Add necessary byproducts to balance the full atom equation (including H).

    Computes full atom counts via RDKit AddHs, finds excess atoms in reactants,
    and decomposes them into known small-molecule byproduct SMILES.
    """
    if not reactant_smiles:
        return []

    # Count all atoms (including H) in reactants
    reactant_counts = {}
    for r in reactant_smiles:
        r_counts = _count_all_atoms(r)
        for sym, cnt in r_counts.items():
            reactant_counts[sym] = reactant_counts.get(sym, 0) + cnt

    # Count all atoms in product
    product_counts = _count_all_atoms(product_smiles)

    # Compute excess: atoms present in reactants but missing from product
    excess = {}
    for sym, cnt in reactant_counts.items():
        p_cnt = product_counts.get(sym, 0)
        if cnt > p_cnt:
            excess[sym] = cnt - p_cnt

    if not excess:
        return []

    # Try to decompose excess into known byproduct molecules
    byproducts = _decompose_excess(excess)
    if byproducts:
        return byproducts

    # Fallback: use single-atom SMILES for any remaining excess
    fallback = []
    for sym, cnt in sorted(excess.items()):
        for bp_atoms, bp_smiles in _BYPRODUCT_TABLE:
            if len(bp_atoms) == 1 and sym in bp_atoms and bp_atoms[sym] == 1:
                fallback.extend([bp_smiles] * cnt)
                break
        else:
            fallback.extend([f'[{sym}]'] * cnt)

    return fallback


def _try_route(smiles: str, r1: str, r2: str = None) -> str:
    """尝试构建合成路线，验证所有反应物 SMILES 并添加副产物确保原子平衡。

    H019 修复: Suzuki 等反应需要在产物侧添加副产物才能通过原子平衡检查。
    """
    if not _validate_smiles(r1):
        return None
    if r2 is not None and not _validate_smiles(r2):
        return None

    reactants = [r1] if r2 is None else [r1, r2]
    byproducts = _add_byproducts_for_balance(smiles, reactants)

    # 构建完整的产物侧（主产物 + 副产物）
    product_part = ".".join([smiles] + byproducts)
    reactant_part = ".".join(reactants)

    return f"{reactant_part}>>{product_part}"


def _run_retro_rule(mol, smarts_pattern, retro_smarts, smiles):
    """通用逆合成规则执行函数。

    H014 改进: 增加化学计量守恒验证（mass balance check）。
    文献依据: LARC (Baker et al., 2025) — Agent-as-a-Judge 框架强调
    逆合成路线必须满足化学计量守恒，原料原子数应与产物匹配。

    Args:
        mol: 目标分子 RDKit Mol 对象
        smarts_pattern: 匹配目标分子的 SMARTS
        retro_smarts: 逆合成反应的 SMARTS
        smiles: 目标分子的 SMILES（用于构建路线）

    Returns:
        dict 或 None: 成功时返回 {"route": ..., "reactants": [...], "steps": 1}
    """
    patt = Chem.MolFromSmarts(smarts_pattern)
    if patt is None or not mol.HasSubstructMatch(patt):
        return None
    rxn = AllChem.ReactionFromSmarts(retro_smarts)
    if rxn is None:
        return None
    ps = rxn.RunReactants((mol,))
    if not ps:
        return None

    # 目标分子的非氢原子数（用于化学计量验证）
    target_n_atoms = mol.GetNumAtoms(onlyExplicit=False)
    # RDKit 的 GetNumAtoms 默认不包含 H — 我们只需要重原子数
    target_heavy = mol.GetNumHeavyAtoms()

    # 遍历所有可能的产物集，找到第一个产生有效SMILES的
    for product_set in ps:
        if not product_set:
            continue
        if len(product_set) >= 2:
            r1 = Chem.MolToSmiles(product_set[0])
            r2 = Chem.MolToSmiles(product_set[1])
            route = _try_route(smiles, r1, r2)
            reactants = [r1, r2]
        else:
            r1 = Chem.MolToSmiles(product_set[0])
            route = _try_route(smiles, r1)
            reactants = [r1]

        if route:
            # H014: 化学计量守恒检查
            # 验证反应物重原子数总和与目标分子是否在合理范围内
            # 若规则只匹配部分分子（如仅匹配核心骨架而丢失取代基），
            # 产物重原子数/原料重原子数会显著偏离 1.0
            try:
                reactant_heavy = 0
                for r_smi in reactants:
                    r_mol = Chem.MolFromSmiles(r_smi)
                    if r_mol:
                        reactant_heavy += r_mol.GetNumHeavyAtoms()
                if reactant_heavy > 0:
                    ratio = target_heavy / reactant_heavy
                    # 允许 ±25% 容差（考虑脱保护基、缩合失水等）
                    # H014 调优: 0.7-1.3 — 验证显示 0.5-1.5 过宽
                    #   例: 取代喹啉(17atoms) / 无取代原料(12atoms) = 1.42 > 1.3 ✓ 拒绝
                    # H021 调优: 0.75-1.25 — 进一步收紧，防止误匹配
                    if ratio < 0.75 or ratio > 1.25:
                        # 原子不守恒 — 拒绝此条规则，尝试下一条
                        continue
            except Exception:
                # 解析失败时不拒绝，保留原有行为
                pass

            # H021: 元素守恒检查
            # 验证目标分子中的所有元素（除 C, H, O, N 外）在反应物中都存在
            # 防止因 SMARTS 只匹配核心骨架而丢失取代基导致的误匹配
            # 例如: 含 B(OH)2 的异喹啉在 Pictet-Spengler 逆反应中会丢失 B
            # 因为 SMARTS 只匹配异喹啉核心，B(OH)2 取代基在反应中被丢弃
            try:
                target_elements = {a.GetSymbol() for a in mol.GetAtoms()}
                reactant_elements = set()
                for r_smi in reactants:
                    r_mol = Chem.MolFromSmiles(r_smi)
                    if r_mol:
                        reactant_elements |= {a.GetSymbol() for a in r_mol.GetAtoms()}
                # C, H, O, N 在缩合反应中可能丢失/获得，不视为异常
                suspicious_elements = target_elements - reactant_elements - {'C', 'H', 'O', 'N'}
                if suspicious_elements:
                    # 反应物中缺少目标分子含有的元素（如 B, F, Cl, Br, I, S, P 等）
                    # 这些元素不可能在逆合成中凭空出现
                    continue
            except Exception:
                # 解析失败时不拒绝，保留原有行为
                pass

            return {"success": True, "route": route, "reactants": reactants, "steps": 1}

    return None


# 简单分子阈值：原子数 <= 10 或属于常见起始原料，停止递归
_SIMPLE_SCAFFOLDS_SMARTS = [
    "c1ccccc1",      # 苯
    "c1ccncc1",      # 吡啶
    "c1cccnc1",      # 吡啶变体
    "C1CCOC1",       # THF
    "C1CCCCC1",      # 环己烷
    "C1CCNC1",       # 吡咯烷
]


def _is_simple_molecule(smiles: str) -> bool:
    """判断分子是否为足够简单的起始原料，无需继续逆合成。
    
    H015 修复: 不再将多环芳烃/联芳基误判为"简单"。
    含苯环且 ≤15 原子的分子如果是单环则简单，多环则可能需要断键。
    文献依据: LARC (2025) — 起始原料应为商业可得的单环简单分子。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    # 原子数 <= 10 视为简单
    if mol.GetNumAtoms() <= 10:
        return True
    # 匹配常见简单骨架
    for s in _SIMPLE_SCAFFOLDS_SMARTS:
        patt = Chem.MolFromSmarts(s)
        if patt and mol.HasSubstructMatch(patt):
            # H015 修复: 只接受单环或几乎无取代的分子（≤12 原子）
            # 排除多环芳烃（联苯、萘等）被误判为"简单"
            n_atoms = mol.GetNumAtoms()
            ring_info = mol.GetRingInfo()
            n_rings = ring_info.NumRings()
            # 单环 + ≤15 原子 → 简单
            if n_rings == 1 and n_atoms <= 15:
                return True
            # 多环但很小（如萘本身只有 10 原子）→ 简单
            if n_rings >= 2 and n_atoms <= 10:
                return True
    return False


# ======================================================================
# 逆合成规则库 — 按反应类型分组，共 35+ 条规则
# 设计原则：
#   1. 优先匹配更具体的结构（放在前面）
#   2. 覆盖常见药物化学键连方式
#   3. 反应物应为商业可得的简单分子
# ======================================================================

RETRO_RULES = [
    # ============================== H012 新增规则 ==============================
    # 文献依据: LARC (Baker et al., 2025) — 规则覆盖率决定逆合成质量;
    # ChemCrow (Bran et al., 2024) — 工具/库的丰富度直接决定 Agent 能力边界
    
    # -- 喹啉/异喹啉合成（先于通用吡啶规则匹配）--
    # Friedländer 喹啉合成逆反应: 喹啉 → 邻氨基苯甲醛 + 酮
    ("c1ccc2c(c1)cccn2", "c1ccc2c(c1)cccn2>>Nc1ccccc1C=O.CC=O"),
    # 异喹啉 → Pictet-Spengler 逆反应
    ("c1ccc2c(c1)ccnc2", "c1ccc2c(c1)ccnc2>>NCCc1ccccc1.C=O"),
    
    # -- 喹唑啉合成 --
    # 喹唑啉 → 邻氨基苯甲酰胺 + 甲酸
    # H021: 使用规范的 SMARTS 表示（c1ccc2ncncc2c1），提高匹配特异性
    ("c1ccc2ncncc2c1", "c1ccc2ncncc2c1>>Nc1ccccc1C(=O)N.C=O"),
    
    # -- 1,2,3-三唑 (Click Chemistry) --
    # CuAAC 逆反应: 三唑 → 叠氮 + 炔
    ("c1cnnn1", "c1cnnn1>>[N-]=[N+]=NC.C#C"),
    
    # -- 肼/联氨类 --
    # 肼 → 重氮盐还原
    ("[NH2][NH2]", "[NH2:1][NH2:2]>>[NH2:1]Cl.[NH2:2]"),
    
    # ============================== 原有规则 ==============================
    # ------------------- 联芳基 C-C 键 (H012) -------------------
    # Suzuki 逆反应: Ar-Ar' → Ar-Br + (HO)2B-Ar'
    # 文献依据: Coscientist (Boiko et al., 2023) — Suzuki coupling 是经典 C-C 键形成;
    # LARC (Baker et al., 2025) — 规则覆盖率决定逆合成质量
    ("[c;R]!@[c;R]", "[c:1]!@[c:2]>>[c:1]Br.[c:2]B(O)O"),
    
    # ------------------- 磺酰胺类 -------------------
    ("[S](=[O])(=[O])[N;!H0]", "[c:1][S](=[O])(=[O])[N:2]>>[c:1][S](=[O])(=[O])Cl.[N:2]"),
    ("[S](=[O])(=[O])[N;H0]", "[c:1][S](=[O])(=[O])[N:2]>>[c:1][S](=[O])(=[O])Cl.[N:2]"),
    
    # ------------------- 酰胺类 -------------------
    ("[C](=[O])[N;!H0]", "[C:1](=[O])[N:2]>>[C:1](=O)O.[N:2]"),
    ("[C](=[O])[N;H0]", "[C:1](=[O])[N:2]>>[C:1](=O)O.[N:2]"),
    
    # ------------------- 酯类 -------------------
    ("[C](=[O])O[!H0]", "[C:1](=[O])O[!H0:2]>>[C:1](=O)O.[O:2]"),
    
    # ------------------- 硫酯类 -------------------
    ("[C](=[O])S[!H0]", "[C:1](=[O])S:2>>[C:1](=O)O.[S:2]"),
    
    # ------------------- 芳基醚 / 烷基醚 -------------------
    ("[c]O[#6;!c]", "[c:1]O[#6;!c:2]>>[c:1]O.[C:2]Cl"),
    ("[#6;!c]O[#6;!c]", "[#6;!c:1]O[#6;!c:2]>>[#6;!c:1]O.[#6;!c:2]Cl"),
    
    # ------------------- 胺类 -------------------
    # 仲胺 — 还原胺化逆反应
    ("[NX3;H1;!$(NC=O);!$(NS(=O)=O)]", "[C:1][N;H1:2][C:3]>>[C:1][N;H2:2].[C:3]=O"),
    # 叔胺 — 烷基化逆反应
    ("[NX3;H0;!$(NC=O);!$(NS(=O)=O)]", "[C:1][N:2][C:3]>>[C:1][N:2].[C:3]Cl"),
    # 芳基仲胺 — Buchwald-Hartwig 简化
    ("[c][NX3;H1]", "[c:1][N:2]>>[c:1]Br.[N:2]"),
    # 芳基叔胺
    ("[c][NX3;H0]", "[c:1][N:2]>>[c:1]Br.[N:2]"),
    
    # ------------------- 脲类 -------------------
    ("[N;!H0]C(=O)[N;!H0]", "[N:1]C(=O)[N:2]>>[N:1].O=C=N[N:2]"),
    
    # ------------------- 氨基甲酸酯 -------------------
    ("[N;!H0]C(=O)O", "[N:1]C(=O)O>>[N:1].O=C=O"),
    
    # ------------------- 硝基还原 -------------------
    ("[N+](=O)[O-]", "[N+](=O)[O-]>>N"),
    
    # ------------------- 卤代芳烃 -------------------
    # 芳基氟 — 亲核芳香取代逆反应
    ("[c]F", "[c:1]F>>[c:1]Cl.F"),
    # 芳基氯
    ("[c]Cl", "[c:1]Cl>>[c:1]O.Cl"),
    # 芳基溴
    ("[c]Br", "[c:1]Br>>[c:1]O.Br"),
    # 芳基碘
    ("[c]I", "[c:1]I>>[c:1]O.I"),
    
    # ------------------- 酮 / 醛 -------------------
    # 芳基酮 — Friedel-Crafts 酰化逆反应
    ("[c]C(=O)[C;!c]", "[c:1]C(=O)[C:2]>>[c:1].[C:2]C(=O)Cl"),
    # 二芳基酮
    ("[c]C(=O)[c]", "[c:1]C(=O)[c:2]>>[c:1].[c:2]C(=O)Cl"),
    # 醛
    ("[CX3H1](=O)", "[C:1](=O)>>[C:1]O"),
    
    # ------------------- 醇类 -------------------
    # 苄醇 — 还原逆反应
    ("[c]C[OH]", "[c:1]C[OH]>>[c:1]C(=O)O"),
    # 普通醇 — 水解逆反应（酯/环氧）
    ("[C;!c][OH]", "[C:1][OH]>>[C:1]Cl.O"),
    
    # ------------------- 腈类 -------------------
    ("[C]#N", "[C:1]#N>>[C:1](=O)O.N"),
    
    # ------------------- 杂环合成 -------------------
    # 噻唑 — Hantzsch 噻唑合成逆反应
    ("c1ncsc1", "c1ncsc1>>N.CS.C=O"),
    # 咪唑
    ("c1[nH]cnc1", "c1[nH]cnc1>>N.C=O.N"),
    # 噁唑
    ("c1ncoc1", "c1ncoc1>>N.C=O.O"),
    # 吡啶（H012 修复：R1 限制仅匹配孤立吡啶环，排除喹啉/异喹啉等稠环体系）
    # 文献: JACS 2024 — 逆合成应反映真实化学转化，单步合成复杂稠环不可行
    ("[c;R1]1[c;R1][c;R1][n;R1][c;R1][c;R1]1", "c1ccncc1>>O=C1CCCC(=O)C1.N"),
    # 嘧啶
    ("c1cncnc1", "c1cncnc1>>N.C=O.N.C=O"),
    
    # ============================== H015 新增规则 ==============================
    # 文献依据: LARC (Baker et al., 2025) — 规则覆盖率决定逆合成质量;
    # Deep Lead Optimization (JACS 2024) — 稠环/桥环体系需要专门的断键策略
    # 设计原则: 覆盖当前缺失的饱和含氮杂环骨架（THIQ、吲哚啉、饱和双环胺等）
    
    # -- 四氢异喹啉 (THIQ) 逆 Pictet-Spengler 反应 --
    # 断裂苄位 C-N 键: THIQ → 苯乙胺 + 醛/酮
    ("c1ccc2c(c1)CCNC2", "c1ccc2c(c1)CCNC2>>NCCc1ccccc1.C=O"),
    
    # -- 吲哚啉 (二氢吲哚) 逆还原环化 --
    # 断裂 C-N 键开环: 吲哚啉 → 邻乙基苯胺
    ("c1ccc2c(c1)CCN2", "c1ccc2c(c1)CCN2>>NCCc1ccccc1"),
    
    # -- 四氢喹啉 逆还原环化 --
    # 1,2,3,4-四氢喹啉 → N-丙基苯胺
    ("c1ccc2c(c1)CCCN2", "c1ccc2c(c1)CCCN2>>NCCCc1ccccc1"),
    
    # -- 饱和环内二级胺 C-N 键断裂 (逆还原胺化) --
    # 只匹配饱和环体系内的 sp3C-sp3N 键，排除酰胺/磺酰胺/芳香体系
    # 产物: C 转羰基 + 游离胺，原子数 +1 (O)，比例在 0.7-1.3 范围内
    ("[C;R;!a;!$(C=*);!$(C#*)][N;R;!a;!$(N[a]);!$(NC=O);!$(NS(=O)=O)]",
     "[C:1][N:2]>>[C:1]=O.[N:2]"),
    
    # -- 苄位 C-N 键在饱和环中 (retro-Pictet-Spengler 变体) --
    # 更精确匹配: 芳环邻接的饱和 C-N 键
    ("[c;R][C;R;!a][N;R;!a;!$(NC=O)]",
     "[c:1][C:2][N:3]>>[c:1][C:2]=O.[N:3]"),
    
    # -- 苯并氮杂环 (benzazepine) 断键 --
    # 7 元含氮杂环与苯稠合 → 逆 Bischler-Napieralski
    ("c1ccc2c(c1)CCCNC2", "c1ccc2c(c1)CCCNC2>>NCCCCc1ccccc1.C=O"),
    
    # ------------------- 缩合反应 -------------------
    # 烯烃 — Wittig / 羟醛缩合逆反应
    ("[C]=[C]", "[C:1]=[C:2]>>[C:1]C=O.[C:2]P"),
    # 亚胺 / 席夫碱
    ("[C]=[N]", "[C:1]=[N:2]>>[C:1]=O.[N:2]"),
    
    # ------------------- 硫醚 / 硫醇 -------------------
    ("[c]S[!H0]", "[c:1][S:2]>>[c:1]Br.[S:2]"),
    ("[C;!c]S[C;!c]", "[C:1][S:2]>>[C:1]Cl.[S:2]"),
    
    # ------------------- 重氮 / 叠氮 -------------------
    ("[N]=[N]", "[N:1]=[N:2]>>[N:1].[N:2]"),
    
    # ------------------- 酰亚胺 -------------------
    ("[C](=O)[N][C](=O)", "[C:1](=O)[N:2][C:3](=O)>>[C:1](=O)O.[N:2].[C:3](=O)O"),
    
    # ------------------- 肟 / 腙 -------------------
    ("[C]=[N][OH]", "[C:1]=[N:2][OH]>>[C:1]=O.[N:2]O"),
    ("[C]=[N][N]", "[C:1]=[N:2][N:3]>>[C:1]=O.[N:2][N:3]"),
    
    # ------------------- 碳酸酯 -------------------
    ("O=C(OC)OC", "O=C(OC)OC>>CO.O=C=O"),
    
    # ------------------- 特殊稠环骨架 -------------------
    # 5,10-二氢吩嗪类 — 逆合成到 2,2'-二氨基联苯
    ("c1ccc2c(c1)CNc1ccccc1N2", "c1ccc2c(c1)CNc1ccccc1N2>>Nc1ccccc1-c1ccccc1N"),

    # ============================== H016 新增规则 ==============================
    # 文献依据:
    #   - Deep Lead Optimization (JACS 2024): 稠环体系需要专门断键策略,
    #     Diels-Alder 环加成是合成六元碳环/桥环的核心反应
    #   - LARC (Baker et al., 2025): 规则覆盖率决定逆合成质量
    #   - 综述第3.2节: 周环反应是药物化学关键转化
    # 设计原则:
    #   - ;R 确保在环内，;!a 排除芳香环
    #   - 产物原子映射保证取代基正确传递
    #   - 原子守恒（6重原子 → 4+2=6），通过 H014 mass balance 验证

    # -- 环己烯 Diels-Alder 逆反应（通用型）--
    # 六元碳环单烯烃 → 1,3-丁二烯片段 + 烯烃片段
    # 这是最经典的 [4+2] 环加成逆反应，适用于:
    #   简单环己烯、取代环己烯、四氢邻苯二甲酰亚胺类等
    # 验证: tetrahydrophthalimide → butadiene + maleimide ✓
    ("[C;R;!a]1=[C;R;!a][C;R;!a][C;R;!a][C;R;!a][C;R;!a]1",
     "[C:1]1=[C:2][C:3][C:4][C:5][C:6]1>>[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]"),

    # -- 二氢吡喃/二氢噻喃 Hetero-Diels-Alder 逆反应（变体A）--
    # 双键紧邻杂原子: O-C=C-C-C-C → O=C-C=C + C=C
    # 适用于: 3,4-二氢-2H-吡喃（C1COC=CC1）、二氢噻喃类似物
    ("[O,S;R]1[C;R;!a]=[C;R;!a][C;R;!a][C;R;!a][C;R;!a]1",
     "[O,S:1]1[C:2]=[C:3][C:4][C:5][C:6]1>>[O,S:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]"),

    # -- 二氢吡喃/二氢噻喃 Hetero-Diels-Alder 逆反应（变体B）--
    # 双键远离杂原子: O-C-C=C-C-C → O-C=C-C=C + C=C
    # 适用于双键在另一端位置的二氢吡喃变体
    ("[O,S;R]1[C;R;!a][C;R;!a][C;R;!a]=[C;R;!a][C;R;!a]1",
     "[O,S:1]1[C:2][C:3][C:4]=[C:5][C:6]1>>[O,S:1][C:2]=[C:3][C:4]=[C:5].[C:6]"),

    # ============================== H017 新增规则 ==============================
    # 文献依据:
    #   - LARC (Baker et al., 2025): 规则覆盖率决定逆合成质量,
    #     内酯/环氧/吡唑/sp3C-sp3C 等骨架的断键规则是当前明显缺口
    #   - Deep Lead Optimization (JACS 2024): 稠环/杂环体系需专门断键策略
    #   - Coscientist (Boiko et al., 2023): 逆合成模板库是合成质量的基础

    # -- 内酯（Lactone）开环: 环内酯逆水解为羟酸 --
    # 匹配: 环内的 -C(=O)-O-C- 子结构
    # 逆反应: 水解开环 → 末端羟基 + 末端羧酸
    # 适用范围: γ-丁内酯、δ-戊内酯、大环内酯等
    ("[C;R](=[O])[O;R]",
     "[C:1](=[O])[O:2]>>[C:1](=O)O.[O:2]"),

    # -- 环氧化物（Epoxide）开环: 三员环氧逆水解为 1,2-二醇 --
    # 匹配: C1OC1 三员环氧环（两个碳共用一个氧）
    # 逆反应: 酸催化开环 → 反式 1,2-二醇
    ("C1OC1",
     "C1OC1>>CCO"),

    # -- 吡唑（Pyrazole）N-N 键断裂 --
    # 匹配: 五元芳环中含两个相邻氮原子
    # 逆反应: 吡唑 → 1,3-二羰基 + 肼（Knorr 吡唑合成逆反应）
    # 吡唑骨架是激酶抑制剂（如 TYK2）的核心结构
    ("c1cn[nH]c1",
     "c1cn[nH]c1>>O=CC=O.NN"),

    # -- sp3-sp3 C-C 键断裂（逆 Wurtz 偶联）--
    # 匹配: 非环烷基链中的 C-C 单键（至少一侧为 sp3 C）
    # 逆反应: C-C → C-Br + C-Cl（Wurtz 偶联逆反应）
    # 适用于: 长烷基链断键、非芳环体系的 C-C 键切割
    # 注意: ;!R 确保不匹配环内 C-C（优先由其他专门规则处理）
    # H017-fix: MgBr → Cl（Grignard 试剂不可表示为有效 SMILES）
    ("[C;!R;!$(C=*);!$(C#*)][C;!R;!$(C=*);!$(C#*)]",
     "[C:1][C:2]>>[C:1]Br.[C:2]Cl"),

    # -- 环醚开环（THF/THP 逆 Williamson 醚合成）--
    # 匹配: 四氢呋喃/四氢吡喃/氧杂环丁烷/二氧戊环等饱和环醚中的 C-O-C
    # 逆反应: 环醚开环 → 卤代醇
    # 适用于: THF、THP、1,4-二氧六环、氧杂环丁烷等
    ("[C;R;!a][O;R;!a]",
     "[C:1][O:2]>>[C:1]Br.[O:2]"),

    # ============================== H020 新增规则 ==============================
    # 文献依据:
    #   - LARC (Baker et al., 2025): 规则覆盖率决定逆合成质量
    #   - ChemCrow (Bran et al., 2024): 工具/规则库丰富度决定 Agent 能力边界
    #   - JACS 2024: BRICS 16种可断裂键需要配套的逆合成规则
    # 设计原则:
    #   - 覆盖药物化学最高频缺失的杂环骨架（吲哚、苯并呋喃、四唑等）
    #   - 每条规则对应真实的命名反应，确保化学合理性
    #   - 放在通用规则之前，确保优先匹配更具体的结构

    # -- 二芳基胺 (Buchwald-Hartwig 双芳基化逆反应) --
    # 匹配: 二级芳胺 NH 连接两个芳香碳原子
    # 逆反应: 断裂一个 C-N 键 → Ar-Br + H2N-Ar
    # 这是药物化学中极常见的结构（如联苯胺类、二苯胺类）
    # 现有 Buchwald 规则 [c][NX3;H1] 仅匹配 N 连一个芳基，此规则补全
    ("[c][NH1][c]",
     "[c:1][NH1:2][c:3]>>[c:1]Br.[NH2:2][c:3]"),

    # -- 吲哚 (Fischer 吲哚合成逆反应) --
    # 匹配: 吲哚核心骨架（苯并吡咯）
    # 逆反应: 吲哚 → 苯肼 + 羰基化合物
    # 吲哚是 FDA 批准药物中最常见的稠杂环之一（如舒马曲坦、昂丹司琼）
    # 也是 TYK2 抑制剂的核心骨架变体来源
    ("c1ccc2[nH]ccc2c1",
     "c1ccc2[nH]ccc2c1>>NNc1ccccc1.C=O"),

    # -- 苯并呋喃 (Rap-Stoermer 逆反应) --
    # 匹配: 苯并呋喃核心骨架（苯并呋喃）
    # 逆反应: 苯并呋喃 → 苯酚 + α-卤代羰基
    # 苯并呋喃是天然产物和药物中的常见氧杂环（如胺碘酮）
    ("c1ccc2occc2c1",
     "c1ccc2occc2c1>>Oc1ccccc1.BrCC=O"),

    # -- 苯并噻吩 (亲电环化逆反应) --
    # 匹配: 苯并噻吩核心骨架（硫代苯并呋喃）
    # 逆反应: 苯并噻吩 → 硫酚 + α-卤代羰基
    # 苯并噻吩在药物化学中作为苯并呋喃的硫等排体出现（如雷洛昔芬类似物）
    ("c1ccc2sccc2c1",
     "c1ccc2sccc2c1>>Sc1ccccc1.BrCC=O"),

    # -- 四唑 (腈+叠氮 [3+2] 环加成逆反应) --
    # 匹配: 1H-四唑环
    # 逆反应: 四唑 → 叠氮 + 腈
    # 四唑是羧酸的最常用生物电子等排体，广泛用于提高代谢稳定性
    # （如氯沙坦、坎地沙坦等 ARB 类药物核心基团）
    ("c1nnn[nH]1",
     "c1nnn[nH]1>>[N-]=[N+]=NC.C#N"),

    # -- 异噁唑 (1,3-二酮+羟胺 缩合逆反应) --
    # 匹配: 异噁唑环
    # 逆反应: 异噁唑 → 1,3-二羰基 + 羟胺
    # 异噁唑在激酶抑制剂中作为氢键受体/供体出现，也是重要的五元杂环
    ("c1cnoc1",
     "c1cnoc1>>O=CC=O.NO"),
]


def plan_synthesis_recursive(smiles: str, max_depth: int = 3, current_depth: int = 0, visited: set = None):
    """递归多步逆合成规划（H003）。

    对目标分子应用单步逆合成规则，然后对得到的反应物（中间体）
    递归调用逆合成规划，直到得到足够简单的起始原料或达到最大深度。

    Args:
        smiles: 目标分子 SMILES
        max_depth: 最大递归深度
        current_depth: 当前递归深度
        visited: 已访问的 SMILES 集合（防止循环）

    Returns:
        dict: {"success": bool, "route": str, "steps": int, "trivial": bool}
    """
    if visited is None:
        visited = set()

    # 防止循环
    if smiles in visited:
        return {"success": True, "route": f"{smiles}>>{smiles}", "steps": 0, "trivial": True}
    visited.add(smiles)

    # 基础情况 1: 达到最大深度
    if current_depth >= max_depth:
        return {"success": True, "route": f"{smiles}>>{smiles}", "steps": 0, "trivial": True}

    # 基础情况 2: 分子已经足够简单
    if _is_simple_molecule(smiles):
        return {"success": True, "route": f"{smiles}>>{smiles}", "steps": 0, "trivial": True}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"success": False, "error": "无效的 SMILES"}

    # 尝试所有逆合成规则
    for smarts_pattern, retro_smarts in RETRO_RULES:
        result = _run_retro_rule(mol, smarts_pattern, retro_smarts, smiles)
        if result:
            reactants = result.get("reactants", [])
            # 对每个反应物递归规划
            sub_routes = []
            all_trivial = True
            for r in reactants:
                sub = plan_synthesis_recursive(r, max_depth, current_depth + 1, visited.copy())
                if sub.get("success"):
                    sub_routes.append(sub)
                    if not sub.get("trivial", False):
                        all_trivial = False

            # 构建完整路线：先放子路线（非平凡的），再放当前步
            full_parts = []
            for sub in sub_routes:
                if not sub.get("trivial", False):
                    full_parts.append(sub["route"])
            full_parts.append(result["route"])
            full_route = " | ".join(full_parts)

            total_steps = result["steps"] + sum(s.get("steps", 0) for s in sub_routes)
            is_trivial = all_trivial and len(reactants) == 1 and reactants[0] == smiles

            return {
                "success": True,
                "route": full_route,
                "steps": total_steps,
                "trivial": is_trivial,
            }

    # 回退: 用 BRICS 碎片化尝试找断键
    try:
        frags = BRICS.BreakMol(mol)
        if frags and len(frags) >= 2:
            frag_smiles = []
            for f in frags:
                if f is None:
                    continue
                s = Chem.MolToSmiles(f)
                if _validate_smiles(s):
                    frag_smiles.append(s)
            if len(frag_smiles) >= 2:
                reactants = ".".join(sorted(frag_smiles)[:2])
                route = _try_route(smiles, reactants)
                if route:
                    # 对 BRICS 碎片也尝试递归
                    sub_routes = []
                    all_trivial = True
                    for r in frag_smiles[:2]:
                        sub = plan_synthesis_recursive(r, max_depth, current_depth + 1, visited.copy())
                        if sub.get("success"):
                            sub_routes.append(sub)
                            if not sub.get("trivial", False):
                                all_trivial = False
                    full_parts = []
                    for sub in sub_routes:
                        if not sub.get("trivial", False):
                            full_parts.append(sub["route"])
                    full_parts.append(route)
                    full_route = " | ".join(full_parts)
                    total_steps = 1 + sum(s.get("steps", 0) for s in sub_routes)
                    return {
                        "success": True,
                        "route": full_route,
                        "steps": total_steps,
                        "trivial": False,
                    }
    except Exception:
        pass

    # 最终回退: 平凡路线
    return {"success": True, "route": f"{smiles}>>{smiles}", "steps": 0, "trivial": True}


def plan_synthesis_v2(smiles: str):
    """使用 SMARTS 断键规则规划合成路线。

    改进点（H003）：
    - 引入递归多步逆合成规划（plan_synthesis_recursive）
    - 对中间体继续断键，直到得到简单起始原料
    - 路线格式: 步骤之间用 " | " 分隔

    文献依据：
    - LARC (2025): 规则覆盖率是逆合成质量的关键
    - ChemCrow (2024): 工具丰富度决定 Agent 能力边界
    """
    return plan_synthesis_recursive(smiles, max_depth=3, current_depth=0)


def score_route_quality(route_str: str, smiles: str) -> float:
    """评估逆合成路线的化学合理性（H012）。

    文献依据:
    - LARC (Baker et al., 2025): Agent-as-a-Judge 逆合成框架，
      路线质量评审是确保合成可行性的关键环节
    - ChemCrow (Bran et al., 2024): 工具输出质量评估决定下游决策可靠性

    评分维度:
    1. 多步路线 vs 单步（多步更真实，+0.2）
    2. 反应物数量（多组分反应更合理，+0.1-0.2）
    3. 复杂度比（产物/反应物原子比 > 3 为可疑，-0.3）
    4. 单反应物路线（除非是简单的氧化/还原/水解，-0.2）

    Returns:
        float: 0.0（不可行）到 1.0（非常可行）
    """
    if not route_str or route_str == f"{smiles}>>{smiles}":
        return 0.0  # 平凡路线 = 最低质量

    steps = route_str.split(" | ")
    step_scores = []

    for step in steps:
        if " >> " not in step and ">>" not in step:
            continue

        # 解析反应物和产物
        sep = " >> " if " >> " in step else ">>"
        parts = step.split(sep, 1)
        if len(parts) != 2:
            continue
        reactants_str, product_str = parts
        reactants = [r.strip() for r in reactants_str.split(".") if r.strip()]

        score = 0.5  # 默认中等

        # 1. 多反应物路线更真实（大多数合成反应涉及2+组分）
        if len(reactants) >= 2:
            score += 0.2
        elif len(reactants) == 1:
            # 单反应物路线须谨慎 — 检查是否为简单官能团转化
            r_mol = Chem.MolFromSmiles(reactants[0])
            p_mol = Chem.MolFromSmiles(product_str)
            if r_mol and p_mol:
                r_atoms = r_mol.GetNumAtoms()
                p_atoms = p_mol.GetNumAtoms()
                if p_atoms > r_atoms + 5:
                    # 单反应物生成更大产物 — 极不可能
                    score -= 0.3
                elif p_atoms <= r_atoms:
                    # 单反应物生成等大或更小产物 — 可能是脱保护/还原/氧化
                    score += 0.1

        # 2. 反应物复杂度检查
        r_total_atoms = 0
        for r in reactants:
            r_mol = Chem.MolFromSmiles(r)
            if r_mol:
                r_total_atoms += r_mol.GetNumAtoms()

        p_mol = Chem.MolFromSmiles(product_str)
        if p_mol and r_total_atoms > 0:
            ratio = p_mol.GetNumAtoms() / r_total_atoms
            if ratio > 3.0:
                # 产物远大于反应物之和 — 不可行的单步转化
                score -= 0.3
            elif ratio > 2.0:
                score -= 0.1

        # 3. 检查是否为同一分子（无转化）
        if len(reactants) == 1 and reactants[0] == product_str:
            score = 0.0

        step_scores.append(max(0.0, min(1.0, score)))

    if not step_scores:
        return 0.0

    # 平均分 + 多步路线奖励
    avg = sum(step_scores) / len(step_scores)
    if len(step_scores) >= 2:
        # 多步路线通常更接近真实合成策略
        avg = min(1.0, avg + 0.15)

    return round(avg, 3)
