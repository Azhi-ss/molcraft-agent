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
from src import config as _config
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
        "run_id": os.environ.get("MOLCRAFT_RUN_ID", "unknown"),
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
    generator: str = Field(
        default="mutate",
        description="生成器选择: mutate(RDKit变异,默认), diffusion(PocketXMol扩散模型,需GPU), hybrid(扩散+RDKit混合)",
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
                generator=params.generator,
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
                    "generator": params.generator,
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
                "run_id": os.environ.get("MOLCRAFT_RUN_ID", "unknown"),
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


def _detect_binding_pocket(residues: list) -> list:
    """无配体时的口袋中心检测。

    策略：计算 N 端 1/3 残基的质心。
    对大多数药物靶点（激酶、GPCR、核受体等），配体结合口袋位于
    蛋白质的 N 端结构域或 N/C 域界面处的 N 端侧，而非两域几何中点。

    这是对 N/C-lobe 中点算法的修正——后者在激酶上偏差 10-14 Å。
    """
    # 收集 CA 原子坐标
    ca_coords = []
    for r in residues:
        ca = None
        for a in r:
            if a.get_name() == 'CA':
                ca = a.get_coord()
                break
        if ca is not None:
            ca_coords.append(ca)

    if not ca_coords:
        return [0.0, 0.0, 0.0]

    coords = np.array(ca_coords)
    n_total = len(coords)

    # N 端 1/3 残基的质心（覆盖大多数蛋白的配体结合域）
    n_third = max(n_total // 3, 10)
    n_term_coords = coords[:n_third]
    n_center = n_term_coords.mean(axis=0)

    # C 端 2/3 残基质心
    c_term_coords = coords[n_third:]
    c_center = c_term_coords.mean(axis=0) if len(c_term_coords) > 0 else n_center

    # 口袋位于 N 端域和 C 端域的交界处，但偏 N 端侧
    # 对于激酶：ATP 口袋在 N-lobe 内，距 N-lobe 质心 ~3 Å
    # 对于 GPCR：正构位点在 N 端跨膜螺旋束内
    # 取 N 端质心 + 向 C 端偏移 ~1/4 的距离
    domain_vec = c_center - n_center
    pocket_center = n_center + domain_vec * 0.25

    return [round(float(c), 2) for c in pocket_center]


def _estimate_box_size(residues: list, center: list, ligand_atoms: list = None) -> list:
    """根据口袋/配体几何自动估算对接盒子大小。

    策略：
    1. 有共晶配体：配体包围盒 + 6 Å margin，下限 18 Å，上限 35 Å
    2. 无配体（apo）：N 端 1/3 残基 CA 散布范围 + 8 Å margin
    3. 默认：25 Å（大多数药物靶点的合理值）
    """
    min_size, max_size = 18.0, 35.0
    default_size = [25.0, 25.0, 25.0]
    center_np = np.array(center)

    if ligand_atoms:
        # 有配体：以配体包围盒为基础
        lig_coords = np.array(ligand_atoms)
        extents = lig_coords.max(axis=0) - lig_coords.min(axis=0)
        margin = 6.0
        size = [round(float(max(s + margin * 2, min_size)), 1) for s in extents]
        # 各维度的上限
        size = [min(s, max_size) for s in size]
        return size

    # 无配体：用 N 端 1/3 残基 CA 的空间散布 + margin
    ca_coords = []
    for r in residues:
        for a in r:
            if a.get_name() == 'CA':
                ca_coords.append(a.get_coord())
                break

    if not ca_coords:
        return default_size

    coords = np.array(ca_coords)
    n_third = max(len(coords) // 3, 10)
    n_term = coords[:n_third]

    # N 端域 CA 的散布范围
    spread = n_term.max(axis=0) - n_term.min(axis=0)
    # 口袋通常在 N 端域内部，盒子需要覆盖该域的主要部分
    margin = 8.0
    size = [round(float(max(s + margin * 2, min_size)), 1) for s in spread]
    size = [min(s, max_size) for s in size]
    return size


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

    # --- 优先用共晶配体重心确定活性位点 ---
    # 常见溶剂/结晶人工产物，跳过
    SOLVENT_NAMES = {
        'HOH', 'WAT', 'EDO', 'GOL', 'DMS', 'SO4', 'PO4', 'ACT', 'IPA',
        'PEG', 'MES', 'HEPES', 'TRS', 'GLC', 'FMT',
        'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'MN', 'FE', 'CU', 'CD', 'CO', 'NI',
    }
    ligand_info = {"found": False, "resname": None, "center": None, "atoms": None}
    for chain_obj in model.get_chains():
        for residue in chain_obj:
            # HETATM 且不是氨基酸、不是溶剂
            if residue.id[0] != ' ' and not is_aa(residue):
                rname = residue.resname.strip().upper()
                if rname in SOLVENT_NAMES:
                    continue
                coords = []
                for atom in residue:
                    if atom.element != 'H':
                        coords.append(atom.get_coord())
                if coords:
                    com = np.array(coords).mean(axis=0)
                    ligand_info = {
                        "found": True,
                        "resname": residue.resname.strip(),
                        "center": [round(float(c), 2) for c in com],
                        "atoms": [[float(x) for x in c] for c in coords],
                    }
                    break
        if ligand_info["found"]:
            break

    if ligand_info["found"]:
        # 用共晶配体重心作为活性位点
        suggested_center = ligand_info["center"]
        center_source = f"共晶配体 {ligand_info['resname']}"
    else:
        # 无配体：用几何口袋检测算法自动定位活性位点
        # 原理: 3D 网格扫描，找到蛋白质表面凹槽最深的位置
        # 适用于任意蛋白靶点，不依赖激酶结构域的先验知识
        suggested_center = _detect_binding_pocket(residues)
        center_source = "几何口袋检测（无共晶配体）"

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
        float(np.linalg.norm(np.array(current_center) - np.array(suggested_center))), 2
    )

    # 自动调整对接坐标
    auto_adjusted = False
    # 自动修正条件:
    # 1. 有共晶配体且偏差 > 2 Å → 自动更新
    # 2. 无配体但偏差 > 10 Å (严重偏移) → 自动更新
    should_auto_fix = (ligand_info["found"] and offset > 2.0) or (not ligand_info["found"] and offset > 10.0)
    if should_auto_fix:
        config_path = Path(__file__).parent.parent / "src" / "config.py"
        try:
            old_line = f"DOCKING_CENTER = {current_center}"
            new_line = f"DOCKING_CENTER = {suggested_center}"
            with open(config_path) as f:
                content = f.read()
            if old_line in content:
                content = content.replace(old_line, new_line)
                with open(config_path, "w") as f:
                    f.write(content)
                # 刷新 import 的 config 缓存
                import importlib
                importlib.reload(config)
                current_center = config.DOCKING_CENTER
                offset = 0.0
                auto_adjusted = True
        except Exception as exc:
            pass

    # 估算推荐盒子大小
    current_size = config.DOCKING_SIZE
    suggested_size = _estimate_box_size(
        residues,
        suggested_center,
        ligand_info.get("atoms") if ligand_info["found"] else None,
    )
    size_offset = round(
        float(np.linalg.norm(np.array(current_size) - np.array(suggested_size))), 2
    )

    # 自动调整盒子大小: 当推荐尺寸与当前值差异 > 10 Å 时自动修正
    size_auto_adjusted = False
    if size_offset > 10.0:
        config_path = Path(__file__).parent.parent / "src" / "config.py"
        try:
            old_line = f"DOCKING_SIZE = {current_size}"
            new_line = f"DOCKING_SIZE = {suggested_size}"
            with open(config_path) as f:
                content = f.read()
            if old_line in content:
                content = content.replace(old_line, new_line)
                with open(config_path, "w") as f:
                    f.write(content)
                importlib.reload(config)
                current_size = config.DOCKING_SIZE
                size_offset = 0.0
                auto_adjusted = True
                size_auto_adjusted = True
        except Exception:
            pass

    size_verdict = "OK" if size_offset < 8.0 else "NEEDS_ADJUSTMENT"

    docking_verdict = "OK" if offset < 8.0 else "NEEDS_ADJUSTMENT"

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
            f"靶点: {protein_info.get('name', 'unknown')} "
            f"({len(seq)} aa), "
            f"活性位点: {center_source}, "
            f"对接坐标偏移 {offset} Å ({docking_verdict}), "
            f"盒子偏移 {size_offset} Å ({size_verdict})"
            + (" [已自动更新 config.py]" if auto_adjusted else "")
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
            "active_site_source": center_source,
            "active_site_cleft": suggested_center,
            "suggested_box_size": suggested_size,
        },
        "docking_check": {
            "current_center": current_center,
            "suggested_center": suggested_center,
            "offset_angstrom": offset,
            "auto_adjusted": auto_adjusted,
            "verdict": docking_verdict,
            "current_size": current_size,
            "suggested_size": suggested_size,
            "size_offset_angstrom": size_offset,
            "size_auto_adjusted": size_auto_adjusted,
            "size_verdict": size_verdict,
        },
        "next_actions": [
            "1. 对接坐标已自动对齐到共晶配体重心" if auto_adjusted else
            f"1. 对接坐标偏差 {offset} Å — 在容差范围内" if docking_verdict == "OK" else
            f"1. ⚠️ 对接坐标偏差 {offset} Å — 建议手动检查",
            f"2. 推荐盒子大小: {suggested_size} (当前: {current_size}, 偏差 {size_offset} Å)" if size_verdict == "NEEDS_ADJUSTMENT" else
            "2. 盒子大小适宜",
            "3. 用蛋白信息指导文献搜索（SearchWeb 搜索蛋白名 + inhibitor/docking）",
            "4. 继续正常迭代流程",
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


# ── Diffusion model (PocketXMol) via HTTP API ──


class DiffusionGenerateParams(BaseModel):
    pdb_path: str = Field(
        default="data/target.pdb",
        description="靶点蛋白 PDB 文件路径",
    )
    n_molecules: int = Field(
        default=20,
        description="生成分子数量，建议 10-50",
    )
    pocket_center: list[float] | None = Field(
        default=None,
        description="口袋中心坐标 [x, y, z]，留空则自动检测",
    )


class DiffusionGenerate(CallableTool2):
    name: str = "diffusion_generate"
    description: str = (
        "用 PocketXMol 扩散模型生成口袋感知分子（需要 GPU 服务器）。"
        "输入蛋白 PDB 和口袋位置，返回 SMILES 列表及药物性质。"
        "输出格式与 generate_molecules 兼容，可直接送入 dock_molecules 或 run_pipeline。"
        "如果 GPU 服务器不可用，会返回错误提示，请改用 generate_molecules。"
    )
    params: type[BaseModel] = DiffusionGenerateParams

    async def __call__(self, params: DiffusionGenerateParams) -> ToolReturnValue:
        api_url = _config.DIFFUSION_API_URL
        if not api_url:
            return ToolError(
                output="",
                message="DIFFUSION_API_URL 未配置。请在 .env 中设置 GPU 服务器地址，或改用 generate_molecules。",
                brief="扩散模型未启用",
            )

        try:
            pdb_path = Path(params.pdb_path)
            if not pdb_path.exists():
                return ToolError(
                    output="",
                    message=f"PDB 文件不存在: {params.pdb_path}",
                    brief="PDB 文件缺失",
                )
            pdb_content = pdb_path.read_text()

            # Build pocket_center from config if not provided
            pocket_center = params.pocket_center
            if pocket_center is None:
                pocket_center = _config.DOCKING_CENTER

            request_body = json.dumps({
                "pdb_content": pdb_content,
                "n_molecules": params.n_molecules,
                "pocket_center": pocket_center,
                "pocket_radius": max(_config.DOCKING_SIZE) / 2.0 if _config.DOCKING_SIZE else 15.0,
            }, ensure_ascii=False).encode("utf-8")

            base_url = api_url.rstrip("/")

            # Phase 1: Submit job
            submit_url = f"{base_url}/generate"
            submit_req = urllib.request.Request(
                submit_url,
                data=request_body,
                headers={"Content-Type": "application/json"},
                method="POST",
            )
            accept_data = await asyncio.to_thread(
                _http_request_json, submit_req, timeout=60
            )
            job_id = accept_data["job_id"]

            # Phase 2: Poll for completion
            poll_interval = 30  # seconds
            timeout = _config.DIFFUSION_TIMEOUT
            elapsed = 0
            while elapsed < timeout:
                await asyncio.sleep(poll_interval)
                elapsed += poll_interval

                status_url = f"{base_url}/job/{job_id}"
                status_req = urllib.request.Request(status_url)
                try:
                    status_data = await asyncio.to_thread(
                        _http_request_json, status_req, timeout=15
                    )
                except urllib.error.URLError:
                    continue  # Transient network issue, keep polling

                job_status = status_data["status"]
                job_elapsed = status_data.get("elapsed_seconds", 0)
                mol_count = status_data.get("molecule_count", 0)

                if job_status == "completed":
                    break
                elif job_status == "failed":
                    return ToolError(
                        output=json.dumps({
                            "status": "error",
                            "error": status_data.get("message", "Unknown"),
                            "hint": f"扩散模型推理失败 (耗时 {job_elapsed:.0f}s)",
                            "fallback": "改用 generate_molecules 生成分子",
                        }),
                        message=status_data.get("message", "Inference failed"),
                        brief="扩散模型推理失败",
                    )

            if elapsed >= timeout:
                return ToolError(
                    output=json.dumps({
                        "status": "error",
                        "error": f"扩散模型推理超时 ({timeout}s)",
                        "hint": "GPU 服务器推理时间超过限制，可能是 PDB 太大或 GPU 过载",
                        "fallback": "改用 generate_molecules 生成分子",
                    }),
                    message=f"Diffusion generation timed out after {timeout}s",
                    brief="扩散模型超时",
                )

            # Phase 3: Fetch result
            result_url = f"{base_url}/job/{job_id}/result"
            result_req = urllib.request.Request(result_url)
            result_data = await asyncio.to_thread(
                _http_request_json, result_req, timeout=30
            )

            molecules = result_data.get("molecules", [])
            gen_time = result_data.get("generation_time_seconds", 0)

            output = {
                "status": "success",
                "summary": f"扩散模型生成了 {len(molecules)} 个口袋感知分子（耗时 {gen_time:.0f}s）",
                "count": len(molecules),
                "molecules": [
                    {
                        "smiles": m["smiles"],
                        "qed": m.get("qed"),
                        "mw": m.get("mw"),
                        "logp": m.get("logp"),
                        "sa_score": m.get("sa_score"),
                        "score": m.get("score"),
                    }
                    for m in molecules
                ],
                "next_actions": [
                    "下一步: 调用 dock_molecules 对这批分子进行对接评估",
                    "或: 直接调用 run_pipeline(generator='hybrid') 混合模式",
                ],
            }
            append_experiment(
                tool="diffusion_generate",
                round_num=get_latest_round() + 1,
                params={"n_molecules": params.n_molecules, "pocket_center": pocket_center},
                result={"output_count": len(molecules), "gen_time": gen_time},
            )
            return ToolOk(output=json.dumps(output, ensure_ascii=False))

        except urllib.error.URLError as exc:
            return ToolError(
                output=json.dumps({
                    "status": "error",
                    "error": str(exc),
                    "hint": "GPU 服务器连接失败。检查 DIFFUSION_API_URL 是否正确，服务器是否运行。",
                    "fallback": "改用 generate_molecules(strategy='mutate') 生成分子",
                }),
                message=str(exc),
                brief="GPU 服务器不可达",
            )
        except Exception as exc:
            return ToolError(
                output=json.dumps({
                    "status": "error",
                    "error": str(exc),
                    "hint": "扩散模型推理失败，可能是 PDB 格式或参数问题",
                    "fallback": "改用 generate_molecules 生成分子",
                }),
                message=str(exc),
                brief="扩散模型推理失败",
            )


def _http_request_json(
    req: urllib.request.Request, timeout: int
) -> dict:
    """Send HTTP request and parse JSON response. Raises on failure."""
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))
