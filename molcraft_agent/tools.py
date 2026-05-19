"""MolCraft Agent 自定义工具定义。

这些工具封装了分子生成、对接、评估和逆合成能力，
供 Kimi CLI Agent 通过 tool_calls 自主调用。

每次工具调用后自动记录到 experiments.jsonl。
"""
import asyncio
import json
import os
import re
import sys
import time
import urllib.request
import urllib.parse
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Blast import NCBIWWW, NCBIXML

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "src"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "tools"))

from pydantic import BaseModel, Field
from kimi_agent_sdk import CallableTool2, ToolError, ToolOk, ToolReturnValue

from generator import generate_molecules
from docking import batch_dock
from synthesis_v2 import plan_synthesis_v2
from evaluator import evaluate_molecule
from receptor import prepare_receptor
from molcraft_agent.experiments import (
    append_experiment,
    get_latest_round,
    get_best_binding_energy,
)
from pipeline import run_evolutionary_pipeline

# ── Stage tracking for begin_stage/end_stage tools ──
_stage_name: str | None = None
_stage_start: float | None = None
_stage_lock = asyncio.Lock()

# ── Molecule persistence for evomap ──

_MOLECULES_JSONL = Path(__file__).resolve().parent.parent / "output" / "molecules.jsonl"


def _persist_molecules(results: list[dict]) -> None:
    """Append pipeline molecule data to molecules.jsonl for evomap."""
    if not results:
        return
    _MOLECULES_JSONL.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "molecules": [
            {
                "smiles": r.get("mol_smiles", ""),
                "be": r.get("binding_energy"),
                "qed": r.get("qed"),
                "trivial": r.get("trivial", False),
            }
            for r in results
        ],
    }
    with open(_MOLECULES_JSONL, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")


class BeginStageParams(BaseModel):
    name: str = Field(description="阶段名称: 诊断, 代码演进, 实验验证, 复盘")


class BeginStage(CallableTool2):
    name: str = "begin_stage"
    description: str = (
        "标记一个阶段的开始。在进入诊断、代码演进、实验验证、复盘四个阶段之一时调用。"
        "必须在对应阶段结束时调用 end_stage 配对。"
    )
    params: type[BaseModel] = BeginStageParams

    async def __call__(self, params: BeginStageParams) -> ToolReturnValue:
        async with _stage_lock:
            global _stage_name, _stage_start
            if _stage_name is not None:
                return ToolError(
                    output="",
                    message=f"阶段「{_stage_name}」尚未结束，不能开始新阶段",
                    brief="阶段嵌套错误",
                )
            _stage_name = params.name
            _stage_start = time.monotonic()
        if hasattr(sys.stdout, "log_stage"):
            sys.stdout.log_stage(name=params.name, status="begun")
        return ToolOk(
            output=json.dumps({
                "stage": params.name,
                "message": f"阶段「{params.name}」开始",
            }, ensure_ascii=False),
        )


class EndStageParams(BaseModel):
    force: bool = Field(
        default=False,
        description="强制结束阶段（用于恢复卡住的状态）。设为 true 时，即使没有正在进行的阶段也不会报错。",
    )


class EndStage(CallableTool2):
    name: str = "end_stage"
    description: str = (
        "标记当前阶段的结束。必须在 begin_stage 之后调用，自动计算耗时并记录。"
        "如果卡住了（begin_stage 后崩溃了），可以传 force=true 强制恢复。"
    )
    params: type[BaseModel] = EndStageParams

    async def __call__(self, params: EndStageParams) -> ToolReturnValue:
        async with _stage_lock:
            global _stage_name, _stage_start
            if _stage_name is None:
                if params.force:
                    return ToolOk(
                        output=json.dumps({
                            "message": "没有正在进行的阶段，force=true 已跳过",
                        }, ensure_ascii=False),
                    )
                return ToolError(
                    output="",
                    message="没有正在进行的阶段，请先调用 begin_stage。如果卡住了，传 force=true",
                    brief="阶段不匹配",
                )
            elapsed = time.monotonic() - _stage_start
            name = _stage_name
            _stage_name = None
            _stage_start = None
        if hasattr(sys.stdout, "log_stage"):
            sys.stdout.log_stage(name=name, status="completed", duration_seconds=round(elapsed, 1))
        return ToolOk(
            output=json.dumps({
                "stage": name,
                "duration_seconds": round(elapsed, 1),
                "message": f"阶段「{name}」完成，耗时 {elapsed:.1f}s",
            }, ensure_ascii=False),
        )


class GenerateParams(BaseModel):
    strategy: str = Field(
        default="mutate",
        description="生成策略: mutate(变异), combine(组合), random(随机)",
    )
    n: int = Field(default=30, description="要生成的分子数量，建议 20-50")
    scaffold: str | None = Field(
        default=None,
        description="可选的种子骨架 SMILES，用于指导生成方向",
    )


class GenerateMolecules(CallableTool2):
    name: str = "generate_molecules"
    description: str = (
        "生成候选药物分子。返回分子 SMILES 列表及其关键性质（QED、MW、LogP）。"
        "运行时间取决于分子数量，通常 5-20 秒。"
        "⚠️ 仅在需要定向探索特定骨架时使用。常规实验请直接用 run_pipeline 工具，"
        "它已包含生成→对接→进化→合成全流程，且默认开启 H002 docking guidance。"
    )
    params: type[BaseModel] = GenerateParams

    async def __call__(self, params: GenerateParams) -> ToolReturnValue:
        try:
            mols = await asyncio.to_thread(
                generate_molecules,
                strategy=params.strategy,
                n_molecules=params.n,
                scaffold=params.scaffold,
            )
            result = {
                "status": "success",
                "summary": f"生成了 {len(mols)} 个候选分子",
                "count": len(mols),
                "molecules": [
                    {
                        "smiles": m["smiles"],
                        "qed": m["qed"],
                        "mw": m["mw"],
                        "logp": m["logp"],
                        "sa_score": m["sa_score"],
                    }
                    for m in mols
                ],
                "next_actions": [
                    "下一步: 调用 dock_molecules 对这批分子进行对接评估",
                    "或: 直接调用 run_pipeline 替代手动流程（推荐，它自带 docking guidance）",
                ],
            }
            append_experiment(
                tool="generate_molecules",
                round_num=get_latest_round() + 1,
                params={"strategy": params.strategy, "n": params.n, "scaffold": params.scaffold},
                result={"output_count": len(mols), "top_qed": mols[0]["qed"] if mols else None},
            )
            return ToolOk(output=json.dumps(result, ensure_ascii=False))
        except Exception as exc:
            return ToolError(
                output=json.dumps({
                    "status": "error",
                    "error": str(exc),
                    "hint": "RDKit 错误通常意味着无效 SMILES——检查种子骨架是否合法",
                    "retry": "用更简单的 scaffold 或 strategy=random 重试",
                }),
                message=str(exc),
                brief="分子生成失败",
            )


class DockParams(BaseModel):
    smiles_list: list[str] = Field(
        description="要对接的分子 SMILES 列表，建议一次不超过 25 个",
    )


class DockMolecules(CallableTool2):
    name: str = "dock_molecules"
    description: str = (
        "批量对分子进行分子对接，计算与靶点的结合自由能（binding_energy，kcal/mol）。"
        "结合能越低（越负）越好，<-7 为优秀。运行时间与分子数量成正比，每个分子约 5-15 秒。"
        "⚠️ 仅在需要评估特定分子时使用。常规实验请直接用 run_pipeline，"
        "它自动完成对接且内置 H002 docking guidance（已验证 +0.4~0.8 kcal/mol）。"
    )
    params: type[BaseModel] = DockParams

    async def __call__(self, params: DockParams) -> ToolReturnValue:
        try:
            await asyncio.to_thread(prepare_receptor)
            mols = [{"smiles": s} for s in params.smiles_list]
            results = await asyncio.to_thread(batch_dock, mols)
            successful = [r for r in results if r.get("success")]
            successful.sort(key=lambda x: x.get("binding_energy", 999))
            top_be = [r.get("binding_energy") for r in successful[:5]]
            output = {
                "total_submitted": len(params.smiles_list),
                "successful": len(successful),
                "failed": len(params.smiles_list) - len(successful),
                "best_be": top_be[0] if top_be else None,
                "top_results": [
                    {
                        "smiles": r["smiles"],
                        "binding_energy": r.get("binding_energy"),
                        "qed": r.get("qed"),
                    }
                    for r in successful[:10]
                ],
            }
            # 自动记录实验
            append_experiment(
                tool="dock_molecules",
                round_num=get_latest_round(),
                params={"n": len(params.smiles_list)},
                result={
                    "successful": len(successful),
                    "failed": len(params.smiles_list) - len(successful),
                    "best_be": top_be[0] if top_be else None,
                    "top5_be": top_be,
                },
            )
            return ToolOk(output=json.dumps(output, ensure_ascii=False))
        except Exception as exc:
            return ToolError(
                output=json.dumps({
                    "status": "error",
                    "error": str(exc),
                    "hint": "对接失败可能是受体文件损坏或 SMILES 无法转换——先检查 prepare_receptor 是否成功",
                    "retry": "重新运行 prepare_receptor 然后重试",
                }),
                message=str(exc),
                brief="分子对接失败",
            )


class SynthesizeParams(BaseModel):
    smiles: str = Field(description="目标分子的 SMILES 字符串")


class PlanSynthesis(CallableTool2):
    name: str = "plan_synthesis"
    description: str = (
        "为目标分子规划逆合成路线，返回 SMILES>>SMILES 格式的反应路线。"
        "支持酰胺、磺酰胺、酯、醚等常见反应类型。"
    )
    params: type[BaseModel] = SynthesizeParams

    async def __call__(self, params: SynthesizeParams) -> ToolReturnValue:
        try:
            result = await asyncio.to_thread(plan_synthesis_v2, params.smiles)
            output = {
                "smiles": params.smiles,
                "success": result.get("success"),
                "route": result.get("route"),
                "steps": result.get("steps"),
            }
            # 自动记录实验
            append_experiment(
                tool="plan_synthesis",
                round_num=get_latest_round(),
                params={"smiles": params.smiles},
                result={"route": result.get("route"), "steps": result.get("steps")},
            )
            return ToolOk(output=json.dumps(output, ensure_ascii=False))
        except Exception as exc:
            return ToolError(
                output="",
                message=str(exc),
                brief="逆合成规划失败",
            )


class EvaluateParams(BaseModel):
    smiles: str = Field(description="要评估的分子的 SMILES 字符串")


class EvaluateMolecule(CallableTool2):
    name: str = "evaluate_molecule"
    description: str = (
        "评估单个分子的药物相似性和理化性质，包括 QED、分子量、LogP、"
        "氢键供体/受体数、可旋转键数、SA score 和 Lipinski 五规则通过情况。"
    )
    params: type[BaseModel] = EvaluateParams

    async def __call__(self, params: EvaluateParams) -> ToolReturnValue:
        try:
            result = await asyncio.to_thread(evaluate_molecule, params.smiles)
            output = {
                "smiles": params.smiles,
                "valid": result.get("valid"),
                "qed": result.get("qed"),
                "mw": result.get("mw"),
                "logp": result.get("logp"),
                "tpsa": result.get("tpsa"),
                "sa_score": result.get("sa_score"),
                "lipinski_pass": result.get("lipinski_pass"),
            }
            # 自动记录实验
            append_experiment(
                tool="evaluate_molecule",
                round_num=get_latest_round(),
                params={"smiles": params.smiles},
                result={
                    "qed": result.get("qed"),
                    "mw": result.get("mw"),
                    "logp": result.get("logp"),
                    "sa_score": result.get("sa_score"),
                },
            )
            return ToolOk(output=json.dumps(output, ensure_ascii=False))
        except Exception as exc:
            return ToolError(
                output="",
                message=str(exc),
                brief="分子评估失败",
            )


class RunPipelineParams(BaseModel):
    n_generate: int = Field(
        default=50,
        description="每代生成的分子数量，建议 30-100",
    )
    n_top: int = Field(
        default=10,
        description="最终保留的 top N 分子",
    )
    strategy: str = Field(
        default="mutate",
        description="生成策略: mutate(变异), combine(组合), random(随机)",
    )
    n_generations: int = Field(
        default=2,
        description="进化代数（包括初始代）",
    )
    use_docking_guidance: bool = Field(
        default=True,
        description="⚡ 铁律：必须保持 True。H002 已验证将结合能从 -7.7~-8.1 提升至 -8.56~-8.80 kcal/mol。关闭它来做'对照实验'是错误的——关闭后基线自然退化 0.5 kcal/mol，无法判断新改动的净效应。正确的 A/B 测试：始终开 docking guidance，在其他变量上做对照。",
    )


class RunPipeline(CallableTool2):
    name: str = "run_pipeline"
    description: str = (
        "⚡ 推荐首选：运行完整的进化迭代药物研发流程。"
        "包含 生成分子→对接→进化→逆合成 全流程，默认开启 H002 docking guidance。"
        "运行时间 5-15 分钟。产出的 result.csv 已保存到 output/ 目录。"
        "除非需要定向探索特定骨架，否则不要用 generate_molecules/dock_molecules 逐个调——那些慢且没有 docking guidance。"
        "铁律：use_docking_guidance 必须保持 True（默认），已验证提升结合能 0.4~0.8 kcal/mol。"
    )
    params: type[BaseModel] = RunPipelineParams

    async def __call__(self, params: RunPipelineParams) -> ToolReturnValue:
        try:
            results = await asyncio.to_thread(
                run_evolutionary_pipeline,
                n_generate=params.n_generate,
                n_top=params.n_top,
                strategy=params.strategy,
                n_generations=params.n_generations,
                use_docking_guidance=params.use_docking_guidance,
                output_dir="output",
            )
            energies = [r["binding_energy"] for r in results if r.get("binding_energy") is not None]
            trivial_count = sum(1 for r in results if r.get("trivial"))
            best_be = min(energies) if energies else None
            avg_be = sum(energies) / len(energies) if energies else None

            # Persist molecule data for evomap visualization
            _persist_molecules(results)

            output = {
                "status": "success",
                "summary": f"Pipeline 完成: {len(results)} 分子, 最佳结合能 {best_be} kcal/mol, trivial 比例 {trivial_count}/{len(results)}",
                "molecule_count": len(results),
                "best_binding_energy": best_be,
                "avg_binding_energy": avg_be,
                "trivial_route_count": trivial_count,
                "trivial_ratio": trivial_count / len(results) if results else 0,
                "top_molecules": [
                    {
                        "smiles": r["mol_smiles"],
                        "binding_energy": r.get("binding_energy"),
                        "route": r.get("route", ""),
                    }
                    for r in results
                ],
                "next_actions": [
                    "1. 调用 report_iteration 记录本轮实验（round/hypothesis_id/success/summary）",
                    "2. 将本轮指标与基线 -8.56 kcal/mol 对比",
                    "3. 如果结合能提升：判定 ACCEPTED，记录改动；如果下降：判定 REJECTED，回退代码",
                    "4. 检查 trivial 比例是否 > 30%，如果是，下一轮考虑扩充逆合成规则库",
                ],
            }
            append_experiment(
                tool="run_pipeline",
                round_num=get_latest_round() + 1,
                params={
                    "n_generate": params.n_generate,
                    "n_top": params.n_top,
                    "strategy": params.strategy,
                    "n_generations": params.n_generations,
                    "docking_guidance": params.use_docking_guidance,
                },
                result={
                    "molecule_count": len(results),
                    "best_be": best_be,
                    "avg_be": avg_be,
                    "trivial_ratio": trivial_count / len(results) if results else 0,
                },
            )
            return ToolOk(output=json.dumps(output, ensure_ascii=False))
        except Exception as exc:
            return ToolError(
                output=json.dumps({
                    "status": "error",
                    "summary": f"Pipeline 运行失败",
                    "error": str(exc),
                    "next_actions": [
                        "1. 检查错误日志中的 traceback",
                        "2. RDKit 错误通常意味着无效 SMILES——检查生成器输出",
                        "3. Vina 错误通常意味着受体文件损坏——运行 prepare_receptor 修复",
                        "4. 网络/API 错误——等待后重试",
                    ],
                }),
                message=str(exc),
                brief="进化管道运行失败",
            )


class ReportIterationParams(BaseModel):
    round_num: int = Field(description="当前迭代轮次（从1开始）")
    hypothesis_id: str = Field(description="本轮验证的假设ID，如 H001")
    success: bool = Field(description="假设是否被验证成功")
    summary: str = Field(
        description="迭代结果摘要，包含关键指标变化（如结合能、QED、trivial route比例）",
    )


class ReportIteration(CallableTool2):
    name: str = "report_iteration"
    description: str = (
        "报告当前迭代已完成。每完成一轮「文献→诊断→改代码→实验验证」的完整闭环后必须调用一次，"
        "main.py 通过此工具的调用来统计迭代次数并控制运行终止。"
    )
    params: type[BaseModel] = ReportIterationParams

    async def __call__(self, params: ReportIterationParams) -> ToolReturnValue:
        try:
            record = {
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "round": params.round_num,
                "hypothesis_id": params.hypothesis_id,
                "success": params.success,
                "summary": params.summary,
            }
            log_path = Path("docs/iteration_log.jsonl")
            log_path.parent.mkdir(exist_ok=True)
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(json.dumps(record, ensure_ascii=False) + "\n")

            output = {
                "status": "recorded",
                "round": params.round_num,
                "message": f"第 {params.round_num} 轮迭代已记录",
            }
            # 假设验证完成时记录 hypothesis_validation 事件
            if hasattr(sys.stdout, "log_hypothesis_validation"):
                sys.stdout.log_hypothesis_validation(
                    hypothesis_id=params.hypothesis_id,
                    success=params.success,
                    conclusion=params.summary,
                )
            return ToolOk(output=json.dumps(output, ensure_ascii=False))
        except Exception as exc:
            return ToolError(
                output="",
                message=str(exc),
                brief="迭代记录失败",
            )


AA_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


class IdentifyTargetParams(BaseModel):
    pdb_path: str = Field(
        default="data/target.pdb",
        description="靶点蛋白 PDB 文件路径",
    )


class IdentifyTarget(CallableTool2):
    name: str = "identify_target"
    description: str = (
        "分析靶点蛋白：提取序列、识别蛋白身份（UniProt/BLAST）、定位活性位点、"
        "验证对接盒子坐标是否对准活性口袋。这是阶段一「文献解析」的第一步，"
        "确保所有后续对接实验基于正确的活性位点坐标。运行时间 5-30 秒。"
    )
    params: type[BaseModel] = IdentifyTargetParams

    async def __call__(self, params: IdentifyTargetParams) -> ToolReturnValue:
        try:
            result = await asyncio.to_thread(_identify_target_impl, params.pdb_path)
            return ToolOk(output=json.dumps(result, ensure_ascii=False))
        except Exception as exc:
            return ToolError(
                output=json.dumps({
                    "status": "error",
                    "error": str(exc),
                    "hint": "PDB 解析失败——检查 target.pdb 是否存在且格式正确",
                    "retry": "确认 data/target.pdb 路径无误后重试",
                }),
                message=str(exc),
                brief="靶点蛋白分析失败",
            )


def _identify_target_impl(pdb_path: str) -> dict:
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB 文件不存在: {pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", str(pdb_path))
    model = structure[0]

    chains = list(model.get_chains())
    chain = chains[0]

    residues = [r for r in chain if is_aa(r)]
    seq = "".join(AA_CODES.get(r.resname, "X") for r in residues)

    # Structural analysis: N-lobe vs C-lobe, hinge region
    heavy_atoms = []
    for r in residues:
        for a in r:
            if a.element != "H":
                heavy_atoms.append(a.get_coord())
    all_coords = np.array(heavy_atoms)
    protein_center = [float(x) for x in all_coords.mean(axis=0)]

    n_lobe_res = residues[:85] if len(residues) >= 85 else residues[:len(residues)//2]
    n_atoms = []
    for r in n_lobe_res:
        for a in r:
            if a.element != "H":
                n_atoms.append(a.get_coord())
    n_center = [float(x) for x in np.array(n_atoms).mean(axis=0)] if n_atoms else protein_center

    c_lobe_res = residues[85:] if len(residues) >= 85 else residues[len(residues)//2:]
    c_atoms = []
    for r in c_lobe_res:
        for a in r:
            if a.element != "H":
                c_atoms.append(a.get_coord())
    c_center = [float(x) for x in np.array(c_atoms).mean(axis=0)] if c_atoms else protein_center

    cleft_center = [float(x) for x in (np.array(n_center) + np.array(c_center)) / 2]
    suggested_center = [round(float(c), 2) for c in cleft_center]

    # UniProt search
    uniprot_result = _search_uniprot(seq)

    # BLAST fallback
    blast_result = None
    if not uniprot_result.get("hit"):
        blast_result = _search_blast(seq)

    # Check current config
    from src import config
    current_center = config.DOCKING_CENTER
    offset = round(
        float(np.linalg.norm(np.array(current_center) - np.array(cleft_center))), 2
    )

    docking_verdict = "OK" if offset < 5.0 else "NEEDS_ADJUSTMENT"

    protein_info = {
        "name": "unknown",
        "species": "unknown",
        "uniprot_id": None,
    }
    if uniprot_result.get("hit"):
        protein_info.update(uniprot_result["hit"])
    elif blast_result:
        protein_info.update(blast_result)

    return {
        "status": "success",
        "summary": (
            f"靶点蛋白: {protein_info.get('name', 'unknown')} "
            f"({len(seq)} aa), 活性位点偏移 {offset} Å, "
            f"对接坐标: {docking_verdict}"
        ),
        "protein": protein_info,
        "sequence": {
            "length": len(seq),
            "sequence": seq,
            "first_60": seq[:60],
            "last_60": seq[-60:],
        },
        "structure": {
            "chains": len(chains),
            "residues": len(residues),
            "protein_center": [round(float(c), 2) for c in protein_center],
            "n_lobe_center": [round(float(c), 2) for c in n_center],
            "c_lobe_center": [round(float(c), 2) for c in c_center],
            "active_site_cleft": suggested_center,
        },
        "docking_check": {
            "current_center": current_center,
            "suggested_center": suggested_center,
            "offset_angstrom": offset,
            "verdict": docking_verdict,
        },
        "next_actions": [
            "1. 如果 docking_check.verdict == 'NEEDS_ADJUSTMENT': 修改 src/config.py 的 DOCKING_CENTER 为 suggested_center",
            "2. 用 identify_target 返回的蛋白信息指导文献搜索（SearchWeb 搜索蛋白名 + inhibitor/docking）",
            "3. 继续阶段一：阅读 papers/ 文献，结合蛋白结构特征提出假设",
        ] if docking_verdict == "NEEDS_ADJUSTMENT" else [
            "1. 对接坐标已验证正确（偏移 < 5 Å），可信任现有对接结果",
            "2. 用蛋白信息指导文献搜索",
            "3. 继续正常迭代流程",
        ],
    }


def _search_uniprot(seq: str) -> dict:
    """Search UniProt REST API for protein identity by sequence fragment."""
    query_seq = seq[:50] if len(seq) > 50 else seq
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f'"{query_seq}"',
        "fields": "accession,protein_name,organism_name,ft_act_site,ft_binding",
        "format": "json",
        "size": "3",
    }
    try:
        req = urllib.request.Request(f"{url}?{urllib.parse.urlencode(params)}")
        req.add_header("User-Agent", "MolCraft-Agent/1.0")
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        results = data.get("results", [])
        if results:
            r = results[0]
            name = r.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "unknown")
            organism = r.get("organism", {}).get("scientificName", "unknown")
            return {
                "hit": {
                    "name": name,
                    "species": organism,
                    "uniprot_id": r.get("primaryAccession", "unknown"),
                }
            }
    except Exception:
        pass
    return {"hit": None}


def _search_blast(seq: str) -> dict | None:
    """NCBI BLAST fallback for protein identification."""
    try:
        handle = NCBIWWW.qblast("blastp", "pdb", seq, hitlist_size=3, expect=0.001)
        records = NCBIXML.parse(handle)
        for rec in records:
            for align in rec.alignments[:1]:
                title = align.title
                for hsp in align.hsps[:1]:
                    identity = round(hsp.identities / hsp.align_length * 100, 1)
                    return {
                        "name": title[:80],
                        "species": "see PDB entry",
                        "identity_pct": identity,
                        "evalue": hsp.expect,
                    }
        handle.close()
    except Exception:
        pass
    return None
