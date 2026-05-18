#!/usr/bin/env python3
"""一键运行完整药物研发流程（进化迭代版）。

改进点（H001）:
- 引入进化式迭代生成：生成 → 对接 → 选择 → 变异 → 再对接
- 文献依据:
  - MOOSE-Chem (Yang et al., 2025): 进化算法导航组合空间
  - Coscientist (Boiko et al., 2023): 基于实验结果的迭代反思
  - MolLEO (Wang et al., 2024b): LLM作为变异和重组算子
"""
import sys
import os
import json
import csv
import tempfile
import uuid
import argparse
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from rdkit import Chem
from generator import generate_molecules, random_mutate_smiles, generate_with_docking_guidance
from docking import batch_dock, dock_molecule, dock_molecule_consensus
from synthesis_v2 import plan_synthesis_v2
from receptor import prepare_receptor
from src.event_schema import MoleculeEvent, MetricsEvent


def log(msg, log_lines):
    ts = datetime.now().isoformat()
    line = f"[{ts}] {msg}"
    log_lines.append(line)
    print(line, file=sys.stderr)


def count_rings(smiles: str) -> int:
    """Count the number of rings in a molecule."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return Chem.rdMolDescriptors.CalcNumRings(mol)



def run_evolutionary_pipeline(
    n_generate=50,
    n_top=10,
    strategy="mutate",
    n_generations=2,
    n_offspring_per_seed=3,
    output_dir="output",
    use_docking_guidance=True,
    logger=None,  # Optional EventLogger instance
):
    """运行进化式迭代药物研发流程。

    Args:
        n_generate: 每代生成的分子数量
        n_top: 最终保留的top N分子
        strategy: 生成策略
        n_generations: 进化代数（包括初始代）
        n_offspring_per_seed: 每个种子产生的变异体数量
        output_dir: 输出目录
        use_docking_guidance: 是否使用 H002 对接引导生成
        logger: 可选的 EventLogger 实例，用于结构化日志输出
    """
    os.makedirs(output_dir, exist_ok=True)
    log_lines = []

    log("=" * 60, log_lines)
    log("MolCraft Agent 进化迭代版开始执行", log_lines)
    log(f"配置: n_generate={n_generate}, n_generations={n_generations}, "
        f"n_offspring={n_offspring_per_seed}, strategy={strategy}, "
        f"docking_guidance={use_docking_guidance}", log_lines)
    log("=" * 60, log_lines)

    # 步骤 1: 准备受体
    log("步骤 1: 准备受体", log_lines)
    prepare_receptor()
    log("受体准备完成", log_lines)

    # 进化迭代
    all_docked = []  # 累积所有对接成功的分子
    current_seeds = None  # 当前代的种子

    for gen in range(n_generations):
        log("-" * 60, log_lines)
        log(f"进化第 {gen + 1}/{n_generations} 代", log_lines)
        log("-" * 60, log_lines)

        if gen == 0:
            # 初始代
            if use_docking_guidance:
                # H002: 使用对接引导生成
                log(f"初始代: 使用对接引导生成 (batch_size=10, n_generations=3)", log_lines)
                mols = generate_with_docking_guidance(
                    docking_fn=dock_molecule,
                    n_molecules=n_generate,
                    batch_size=10,
                    n_generations=3,
                    top_k=5,
                    strategy=strategy,
                )
                log(f"对接引导生成完成，获得 {len(mols)} 个分子", log_lines)
            else:
                log(f"初始代: 生成候选分子 (strategy={strategy}, n={n_generate})", log_lines)
                mols = generate_molecules(strategy=strategy, n_molecules=n_generate)
        else:
            # 后续代：从上一代top种子变异生成
            log(f"第 {gen + 1} 代: 从 top {len(current_seeds)} 种子变异生成", log_lines)
            mols = _generate_offspring(current_seeds, n_offspring_per_seed)

        log(f"本代生成 {len(mols)} 个通过过滤的分子", log_lines)

        # 对接（如果用了 docking_guidance，mols 已经包含对接结果）
        if use_docking_guidance and gen == 0:
            # 初始代已对接，直接使用
            successful = [d for d in mols if d.get("docking_success")]
            successful.sort(key=lambda x: x.get("binding_energy", 999))
            log(f"对接引导生成已包含对接结果，成功 {len(successful)}/{len(mols)}", log_lines)
        else:
            log(f"第 {gen + 1} 代: 分子对接", log_lines)
            docked = batch_dock(mols)
            successful = [d for d in docked if d.get("success")]
            successful.sort(key=lambda x: x.get("binding_energy", 999))
            log(f"对接成功 {len(successful)}/{len(mols)} 个分子", log_lines)

        # 累积到全局池
        all_docked.extend(successful)
        all_docked.sort(key=lambda x: x.get("binding_energy", 999))

        # 选择下一代种子（取本代top，避免过度收敛）
        n_seeds = min(20, len(successful))
        current_seeds = successful[:n_seeds]
        log(f"选择 {n_seeds} 个种子进入下一代", log_lines)
        if successful:
            log(f"本代最佳结合能: {successful[0].get('binding_energy')} kcal/mol", log_lines)
            log(f"全局最佳结合能: {all_docked[0].get('binding_energy')} kcal/mol", log_lines)

    # 最终选择
    log("-" * 60, log_lines)
    log("最终选择: 从所有代中选择 top 分子", log_lines)

    # 去重（按SMILES）
    seen = set()
    unique_docked = []
    for d in all_docked:
        smi = d.get("smiles", "")
        if smi and smi not in seen:
            seen.add(smi)
            unique_docked.append(d)

    unique_docked.sort(key=lambda x: x.get("binding_energy", 999))

    # H009: 共识对接 — 对 top 候选分子用多次独立对接取中位数
    # Coscientist 模式："performing experiments multiple times"
    # 消除单次对接的随机噪声，提升最终排名可靠性
    N_CONSENSUS = min(n_top * 2, len(unique_docked))
    log(f"共识对接: 对 top {N_CONSENSUS} 候选执行 3 次独立对接取中位数...", log_lines)
    consensus_top = unique_docked[:N_CONSENSUS]
    consensus_results = []
    for i, candidate in enumerate(consensus_top):
        smiles = candidate.get("smiles", "")
        if not smiles:
            continue
        cresult = dock_molecule_consensus(smiles)
        if cresult.get("success"):
            candidate["binding_energy"] = cresult["binding_energy"]
            candidate["consensus_std"] = cresult.get("std_energy", 0.0)
            candidate["consensus_n"] = cresult.get("n_success", 0)
            consensus_results.append(candidate)
            log(f"  [{i+1}/{N_CONSENSUS}] {smiles[:40]}... "
                f"median={cresult['binding_energy']:.3f} std={cresult['std_energy']:.3f} "
                f"(n={cresult['n_success']})", log_lines)

    # 按中位数结合能重新排序
    consensus_results.sort(key=lambda x: x.get("binding_energy", 999))

    # 合并回 unique_docked：共识候选用新分数，其余保持不变
    consensus_smiles = {r["smiles"] for r in consensus_results}
    remaining = [d for d in unique_docked if d.get("smiles") not in consensus_smiles]
    unique_docked = consensus_results + remaining

    # H012: 路线质量评分（LARC Agent-as-a-Judge 模式）
    # 对候选池做逆合成规划 → 评分 → 复合排序
    from synthesis_v2 import score_route_quality

    candidate_pool_size = min(n_top * 3, len(unique_docked))
    candidate_pool = unique_docked[:candidate_pool_size]
    log(f"去重后共 {len(unique_docked)} 个独特分子，候选池 {candidate_pool_size} 个进行路线评分", log_lines)

    if not candidate_pool:
        log("错误: 没有分子对接成功", log_lines)
        sys.exit(1)

    # 逆合成规划 + 路线质量评分
    log(f"为候选池 {len(candidate_pool)} 个分子规划合成路线并评分", log_lines)
    scored_candidates = []
    for mol in candidate_pool:
        smiles = mol["smiles"]
        syn = plan_synthesis_v2(smiles)
        route = syn.get("route", f"{smiles}>>{smiles}") if syn.get("success") else f"{smiles}>>{smiles}"
        is_trivial = syn.get("trivial", False) or route == f"{smiles}>>{smiles}"
        quality = score_route_quality(route, smiles)
        scored_candidates.append({
            "mol_smiles": smiles,
            "route": route,
            "binding_energy": mol.get("binding_energy"),
            "qed": mol.get("qed"),
            "trivial": is_trivial,
            "route_quality": quality,
            "syn_steps": syn.get("steps", 0),
        })

        if logger is not None:
            composite = 0.0  # will be recalculated in H012 scoring below
            logger.log_molecule(MoleculeEvent(
                mol_smiles=smiles,
                binding_energy=mol.get("binding_energy") or 0.0,
                composite_score=composite,
                syn_steps=syn.get("steps", 0),
                sa_score=syn.get("sa_score"),
                qed=mol.get("qed") or 0.0,
                rings=count_rings(smiles),
                trivial=is_trivial,
                route_quality=quality,
                docking_std=mol.get("consensus_std"),
            ))

    # H012 复合评分: 0.8×BE_norm + 0.2×route_quality
    energies_all = [c["binding_energy"] for c in scored_candidates if c["binding_energy"] is not None]
    if energies_all:
        e_min, e_max = min(energies_all), max(energies_all)
        if e_max > e_min:
            for c in scored_candidates:
                if c["binding_energy"] is not None:
                    be_norm = (e_max - c["binding_energy"]) / (e_max - e_min)
                else:
                    be_norm = 0.0
                c["composite_score"] = 0.8 * be_norm + 0.2 * c["route_quality"]
        else:
            for c in scored_candidates:
                c["composite_score"] = c["route_quality"]

        # 按复合评分排序（越高越好）
        scored_candidates.sort(key=lambda x: x.get("composite_score", 0), reverse=True)
        log(f"复合评分排序完成 (0.8×BE + 0.2×路线质量)", log_lines)

    # 选择 top N
    final_top = scored_candidates[:n_top]
    results = []
    trivial_count = 0
    for c in final_top:
        if c["trivial"]:
            trivial_count += 1
        results.append({
            "mol_smiles": c["mol_smiles"],
            "route": c["route"],
            "binding_energy": c["binding_energy"],
            "qed": c["qed"],
            "trivial": c["trivial"],
        })
        log(f"  {c['mol_smiles'][:50]}... BE={c['binding_energy']} "
            f"steps={c['syn_steps']} quality={c['route_quality']:.2f} "
            f"composite={c['composite_score']:.3f}"
            f"{' [TRIVIAL]' if c['trivial'] else ''}", log_lines)

    # 统计
    energies = [r["binding_energy"] for r in results if r["binding_energy"] is not None]
    if energies:
        avg_energy = sum(energies) / len(energies)
        log(f"Top {len(results)} 平均结合能: {avg_energy:.3f} kcal/mol", log_lines)
        log(f"最佳结合能: {min(energies):.3f} kcal/mol", log_lines)
        log(f"Trivial route 比例: {trivial_count}/{len(results)} ({trivial_count/len(results)*100:.1f}%)", log_lines)

    # 保存 CSV：先写临时文件，验证后原子替换，避免部分写入文件
    csv_path = os.environ.get("MOLCRAFT_CSV_PATH", os.path.join(output_dir, "result.csv"))
    csv_dir = os.path.dirname(os.path.abspath(csv_path)) or output_dir
    fd, temp_csv = tempfile.mkstemp(prefix="result.", suffix=".csv.tmp", dir=csv_dir)
    os.close(fd)

    try:
        with open(temp_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["mol_smiles", "route"])
            writer.writeheader()
            for row in results:
                writer.writerow({"mol_smiles": row["mol_smiles"], "route": row["route"]})

        # 原子替换：只有完整写入成功后才替换正式文件
        os.replace(temp_csv, csv_path)
        log(f"结果已保存到 {csv_path}", log_lines)
    except Exception as e:
        if os.path.exists(temp_csv):
            os.unlink(temp_csv)
        log(f"CSV 写入失败: {e}", log_lines)
        raise

    # 结构化日志：通过 EventLogger 输出 metrics
    if logger is not None:
        energies = [r["binding_energy"] for r in results if r["binding_energy"] is not None]
        steps = [r.get("syn_steps", 0) for r in results]
        logger.log_metrics(MetricsEvent(
            molecule_count=len(results),
            non_trivial_count=len(results) - trivial_count,
            trivial_count=trivial_count,
            avg_binding_energy=sum(energies) / len(energies) if energies else None,
            min_binding_energy=min(energies) if energies else None,
            avg_syn_steps=sum(steps) / len(steps) if steps else None,
            docking_success_rate=0.0,
        ))

    # 打印摘要
    log("=" * 60, log_lines)
    log("流程完成", log_lines)
    log(f"最佳结合能: {results[0]['binding_energy']} kcal/mol", log_lines)
    log(f"输出文件: {csv_path}", log_lines)
    log("=" * 60, log_lines)

    return results


def _generate_offspring(seeds, n_offspring_per_seed):
    """从种子分子生成变异后代。

    利用对接成功的种子作为起点，通过随机变异产生新分子。
    这模拟了进化算法中的"选择+变异"步骤。
    """
    import random
    from evaluator import evaluate_molecule, passes_filters

    offspring = []
    attempts = 0
    max_attempts = len(seeds) * n_offspring_per_seed * 20

    while len(offspring) < len(seeds) * n_offspring_per_seed and attempts < max_attempts:
        attempts += 1
        seed = random.choice(seeds)
        seed_smiles = seed.get("smiles", "")
        if not seed_smiles:
            continue

        # 变异强度：1-4个突变点
        n_mut = random.randint(1, 4)
        new_smiles = random_mutate_smiles(seed_smiles, n_mut)

        if new_smiles and new_smiles != seed_smiles:
            # 检查是否为单片段（多片段 SMILES 会导致 Meeko 崩溃）
            mol_check = Chem.MolFromSmiles(new_smiles)
            if mol_check is None or len(Chem.GetMolFrags(mol_check)) > 1:
                continue
            props = evaluate_molecule(new_smiles)
            if passes_filters(props):
                offspring.append(props)
    
    return offspring


def main():
    parser = argparse.ArgumentParser(description="运行进化迭代版药物研发流程")
    parser.add_argument("--n-generate", type=int, default=50, help="每代生成分子数量 (默认: 50)")
    parser.add_argument("--n-top", type=int, default=10, help="保留 top N 分子 (默认: 10)")
    parser.add_argument("--strategy", choices=["mutate", "combine", "random"], default="mutate",
                        help="生成策略 (默认: mutate)")
    parser.add_argument("--n-generations", type=int, default=2, help="进化代数 (默认: 2)")
    parser.add_argument("--n-offspring", type=int, default=3, help="每个种子的变异体数量 (默认: 3)")
    parser.add_argument("--output-dir", type=str, default="output", help="输出目录")
    parser.add_argument("--no-docking-guidance", action="store_false", dest="use_docking_guidance",
                        help="禁用 H002 对接引导（不推荐，已验证会退化 0.4~0.8 kcal/mol）")
    args = parser.parse_args()

    run_evolutionary_pipeline(
        n_generate=args.n_generate,
        n_top=args.n_top,
        strategy=args.strategy,
        n_generations=args.n_generations,
        n_offspring_per_seed=args.n_offspring,
        output_dir=args.output_dir,
        use_docking_guidance=args.use_docking_guidance,
    )


if __name__ == "__main__":
    main()
