# 口袋感知3D分子生成模型开源仓库与论文综合调研报告

## 目录
- [1. DiffSBDD及变体/后续模型](#1-diffsbdd及变体后续模型)
- [2. TargetDiff及变体](#2-targetdiff及变体)
- [3. Pocket2Mol及变体](#3-pocket2mol及变体)
- [4. DiffDock分子对接模型](#4-diffdock分子对接模型)
- [5. 2024-2025最新SBDD/分子生成论文和开源仓库](#5-2024-2025最新sbdd分子生成论文和开源仓库)

---

## 1. DiffSBDD及变体/后续模型

### 1.1 DiffSBDD (核心模型)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/arneschneuing/DiffSBDD |
| Stars | 504 |
| 论文标题 | DiffSBDD: Structure-based Drug Design with Equivariant Diffusion Models |
| 发表年份/会议 | Nature Computational Science, 2024 (arXiv: 2210.13695, 2022) |
| 预训练权重 | https://zenodo.org/record/8183747 (8个模型：CrossDocked/Binding MOAD × conditional/joint × Cα/fullatom) |
| 推理脚本 | `generate_ligands.py` (De novo), `inpaint.py` (子结构修复), `optimize.py` (分子优化) |
| 输出格式 | SDF (`.sdf`) |
| SDF→SMILES转换难度 | 低 (RDKit可直接从SDF读取3D构象并转换为SMILES, 一行代码: `Chem.MolFromMolFile('x.sdf')`) |
| 4090 GPU部署需求 | PyTorch+CUDA10.2+PyG, 模型参数量中等(EGNN), 单卡4090完全可运行推理, 训练也可单卡(需调整batch_size) |

**推理示例：**
```bash
python generate_ligands.py checkpoints/crossdocked_fullatom_cond.ckpt \
  --pdbfile example/3rfm.pdb --outfile example/3rfm_mol.sdf \
  --ref_ligand A:330 --n_samples 20
```

**特点：**
- 支持De novo设计、子结构修复(Inpainting)、分子优化三种模式
- 等变扩散模型(Equivariant Diffusion)，基于EGNN
- 有Colab Notebook可直接试用
- 作者提示：可尝试新模型DrugFlow

### 1.2 DiffSBDD-DDPO (RL优化变体)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/X-east/diffsbdd-ddpo |
| Stars | 0 |
| 论文标题 | DiffSBDD-DDPO: Diffusion Model Policy Optimization for Structure-Based Drug Design (推测) |
| 备注 | 非官方fork，仓库为早期版本，基于DiffSBDD原始代码(非更新版) |
| 推理脚本 | 同DiffSBDD的generate_ligands.py |
| 输出格式 | SDF |
| 4090 GPU部署需求 | 同DiffSBDD |

### 1.3 DrugFlow (DiffSBDD作者后续模型, 强烈推荐)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/LPDI-EPFL/DrugFlow |
| Stars | 158 |
| 论文标题 | Multi-domain Distribution Learning for De Novo Drug Design |
| 发表年份/会议 | ICLR 2025 |
| 预训练权重 | https://zenodo.org/records/14919171 (4个模型：DrugFlow, DrugFlow+OOD, FlexFlow, DrugFlow+PA) |
| 推理脚本 | `src/generate.py` |
| 输出格式 | SDF (`--output examples/samples.sdf`) |
| SDF→SMILES转换难度 | 低 |
| 4090 GPU部署需求 | PyTorch Lightning + Flow Matching + Markov Bridge, 单卡4090可推理, 有Docker支持 |

**推理示例：**
```bash
python src/generate.py \
  --protein examples/kras.pdb \
  --ref_ligand examples/kras_ref_ligand.sdf \
  --checkpoint checkpoints/drugflow.ckpt \
  --output examples/samples.sdf
```

**特点：**
- **DiffSBDD官方作者的后续改进模型**
- 结合连续Flow Matching和离散Markov Bridge
- 支持OOD检测(置信度头)、偏好对齐(PA)、侧链柔性(FlexFlow)
- 当前最先进的SBDD生成模型之一

### 1.4 DiffSBDDFactor (非官方变体)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/rwbfd/DiffSBDDFactor |
| Stars | 0 |
| 备注 | 非官方fork, 功能不明 |

---

## 2. TargetDiff及变体

### 2.1 TargetDiff (核心模型)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/guanjq/targetdiff |
| Stars | 339 |
| 论文标题 | 3D Equivariant Diffusion for Target-Aware Molecule Generation and Affinity Prediction |
| 发表年份/会议 | ICLR 2023 |
| 预训练权重 | https://drive.google.com/drive/folders/1-ftaIrTXjWFhw3-0Twkrs5m0yX6CNarz |
| 推理脚本 | `scripts/sample_diffusion.py`, `scripts/sample_for_pocket.py`, `scripts/batch_sample_diffusion.sh` |
| 输出格式 | 输出为meta信息文件(需后续evaluate脚本处理), 最终可通过RDKit转为SDF |
| SDF→SMILES转换难度 | 中 (需先运行evaluate脚本从meta文件提取分子信息, 然后RDKit转换) |
| 4090 GPU部署需求 | PyTorch 1.13+PyG+CUDA11.6, 单卡4090可推理, 训练推荐多卡 |

**推理示例：**
```bash
# 对测试集口袋采样
python scripts/sample_diffusion.py configs/sampling.yml --data_id {i}
# 对自定义PDB口袋采样
python scripts/sample_for_pocket.py configs/sampling.yml --pdb_path examples/xxx.pdb
# 多GPU批量采样
CUDA_VISIBLE_DEVICES=0 bash scripts/batch_sample_diffusion.sh configs/sampling.yml outputs 4 0 0
```

**特点：**
- 同时支持分子生成和亲和力预测(双功能)
- 等变扩散模型，基于EGNN
- 提供完整的Vina Docking评估流水线
- 与Pocket2Mol共享数据集处理流程

### 2.2 DecompDiff (TargetDiff变体, 子结构分解先验)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/bytedance/DecompDiff |
| Stars | 74 |
| 论文标题 | DecompDiff: Diffusion Models with Decomposed Priors for Structure-Based Drug Design |
| 发表年份/会议 | ICML 2023 (字节跳动) |
| 预训练权重 | https://drive.google.com/drive/folders/1JAB5pp25rEM5Wt-i373_rrAyTsLvAACZ (主模型); https://drive.google.com/drive/folders/1QOQOuDxdKkipYygZU9OIQUXqV9C28J5O (beta priors) |
| 推理脚本 | `scripts/sample_diffusion_decomp.py` |
| 输出格式 | meta信息文件 → evaluate脚本处理 → SDF/PDB |
| SDF→SMILES转换难度 | 中 (同TargetDiff流程) |
| 4090 GPU部署需求 | 同TargetDiff环境, 额外需mdtraj+alphaspace2, 单卡4090可推理 |

**推理示例：**
```bash
python scripts/sample_diffusion_decomp.py configs/sampling_drift.yml \
  --outdir $SAMPLE_OUT_DIR -i $DATA_ID --prior_mode ref_prior
```

**特点：**
- 引入子结构分解先验(decomposed priors)显著提升Vina Dock得分
- 支持ref_prior和beta_prior两种模式
- 在Vina Dock(-8.43)和High Affinity(71%)指标上大幅超越TargetDiff和Pocket2Mol

### 2.3 BindDM (TargetDiff系变体)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/YangLing0818/BindDM |
| Stars | 21 |
| 论文标题 | Binding-Adaptive Diffusion Models for Structure-Based Drug Design |
| 发表年份/会议 | AAAI 2024 |
| 预训练权重 | 未明确提供(需自行训练) |
| 推理脚本 | `sample.py` |
| 输出格式 | 同TargetDiff格式(meta文件) |
| SDF→SMILES转换难度 | 中 |
| 4090 GPU部署需求 | 同TargetDiff环境 |

**特点：**
- 利用蛋白质-配体子复合体(subcomplex)增强结合适应性生成
- 数据准备和评估流程与TargetDiff完全兼容

---

## 3. Pocket2Mol及变体

### 3.1 Pocket2Mol (核心模型)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/pengxingang/Pocket2Mol |
| Stars | 401 |
| 论文标题 | Pocket2Mol: Efficient Molecular Sampling Based on 3D Protein Pockets |
| 发表年份/会议 | ICLR 2023 (arXiv: 2205.07249) |
| 预训练权重 | 仓库内ckpt/README.md指引下载 |
| 推理脚本 | `sample.py` (测试集), `sample_for_pdb.py` (自定义PDB), `batch_sample.sh` (批量) |
| 输出格式 | 输出为pickle+内部格式, 可通过评估脚本转为SDF |
| SDF→SMILES转换难度 | 中 (需从内部格式提取3D坐标+原子信息, 再构建RDKit Mol对象) |
| 4090 GPU部署需求 | PyTorch 1.10+PyG 2.0+CUDA 11.3, 单卡4090可推理(推荐taskset限制CPU核心加速) |

**推理示例：**
```bash
# 自定义PDB口袋采样
python sample_for_pdb.py \
  --pdb_path ./example/4yhj.pdb \
  --center " 32.0,28.0,36.0"
# 测试集采样
CUDA_VISIBLE_DEVICES=0 taskset -c 0 python sample.py --data_id 0 --outdir ./outputs
```

**特点：**
- 自回归式采样(逐原子添加), 比扩散模型更高效
- 等变GNN架构
- 支持自定义PDB口袋输入(只需PDB文件+口袋中心坐标)
- 发现单CPU核心采样更快(需taskset -c)

### 3.2 Pocket2Mol-RL (RL微调变体)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/deargen/Pocket2Mol_RL_public |
| Stars | 10 |
| 论文标题 | Fine-tuning Pocket-conditioned 3D Molecule Generation via Reinforcement Learning |
| 发表年份/会议 | ICLR 2024 Workshop? |
| 预训练权重 | 仓库内提供推理数据+模型 |
| 推理脚本 | 内部推理脚本(需pip install当前仓库) |
| 输出格式 | 同Pocket2Mol |
| 4090 GPU部署需求 | 同Pocket2Mol + Docker支持 |

**特点：**
- 通过强化学习(RL)微调Pocket2Mol, 提升生成分子的药化性质
- 不提供RL训练代码, 仅提供推理代码
- 有Docker镜像 `deargen/pocket2mol_rl_public:latest`

### 3.3 Apo2Mol (Pocket2Mol系变体, Apo口袋条件)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/AIDD-LiLab/Apo2Mol |
| Stars | 32 |
| 论文标题 | Apo2Mol: 3D Molecule Generation via Dynamic Pocket-Aware Diffusion Models |
| 发表年份/会议 | AAAI 2026 (arXiv: 2511.14559) |
| 预训练权重 | 仓库内包含 `apo2mol_checkpoint.ckpt` |
| 数据集 | https://huggingface.co/datasets/AIDD-LiLab/Apo2Mol_Dataset |
| 推理脚本 | `sample_split.py` |
| 输出格式 | 内部格式(需eval_split.py评估) |
| SDF→SMILES转换难度 | 中 |
| 4090 GPU部署需求 | PyTorch Lightning + Hydra + W&B, 单卡可推理, 多卡推荐 |

**推理示例：**
```bash
CUDA_VISIBLE_DEVICES=0 python sample_split.py   # 单GPU
CUDA_VISIBLE_DEVICES=0,1,2,3 python sample_split.py  # 多GPU
```

**特点：**
- **首个基于Apo(无配体)蛋白口袋的3D分子生成模型**
- 联合生成配体和Holo口袋构象(解决Apo-Holo口袋形变问题)
- 实用性极强——实际SBDD场景通常只有Apo蛋白结构
- 数据集发布在HuggingFace

---

## 4. DiffDock分子对接模型(评分/对接)

### 4.1 DiffDock (核心对接模型)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/gcorso/DiffDock |
| Stars | 1512 |
| 论文标题 | DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking |
| 发表年份/会议 | ICLR 2023 (arXiv: 2210.01776); DiffDock-L: ICLR 2024 |
| 预训练权重 | 仓库内含, 默认运行DiffDock-L(新版本) |
| 推理脚本 | `python -m inference` + `default_inference_args.yaml` |
| 输出格式 | PDB文件(预测的配体pose) |
| SDF→SMILES转换难度 | 不适用(DiffDock是对接模型, 输入为蛋白PDB+配体SDF/SMILES, 输出为pose PDB, 不需要转SMILES) |
| 4090 GPU部署需求 | conda环境+CUDA, 单卡4090可推理, 有Docker镜像 `rbgcsail/diffdock`, 也可CPU运行(慢) |

**推理示例：**
```bash
python -m inference --config default_inference_args.yaml \
  --protein_ligand_csv data/protein_ligand_example.csv \
  --out_dir results/user_predictions_small
# 单复合物:
python -m inference --protein_path protein.pdb --ligand "COc(cc1)ccc1C#N"
```

**特点：**
- **分子对接评分模型(不是生成模型)**, 输出配体pose和置信度分数
- DiffDock-L(2024)版本显著提升性能和泛化能力
- 支持HuggingFace Spaces在线试用
- 支持PDB文件、ESMFold序列、SMILES、SDF等多种输入
- 有GUI界面 `python app/main.py`
- 输出置信度分数: c>0表示可靠预测

### 4.2 DiffDock-Pocket (口袋级对接+侧链柔性)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/plainerman/DiffDock-Pocket |
| Stars | 38 |
| 论文标题 | DiffDock-Pocket: Diffusion for Pocket-Level Docking with Side Chain Flexibility |
| 预训练权重 | 仓库内含模型权重 |
| 推理脚本 | `inference.py` |
| 输出格式 | PDB文件(预测pose+侧链构象) |
| 4090 GPU部署需求 | 同DiffDock, 单卡4090可推理, 推荐batch_size=20, samples_per_complex=40 |

**推理示例：**
```bash
python inference.py --protein_path example_data/3dpf_protein.pdb \
  --ligand example_data/3dpf_ligand.sdf \
  --batch_size 20 --samples_per_complex 40 \
  --keep_local_structures --save_visualisation
```

**特点：**
- DiffDock的口袋级扩展, 支持侧链柔性
- 适用于计算预测蛋白结构(如ESMFold)
- 输入也可直接用SMILES代替SDF
- 支持能量最小化(`--relax`)

### 4.3 DisCo-DiffDock (离散+连续扩散增强对接)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/gcorso/disco-diffdock |
| Stars | 93 |
| 论文标题 | DisCo-Diff: Enhancing Continuous Diffusion Models with Discrete Latents |
| 发表年份/会议 | ICML 2024 |
| 推理脚本 | `python -m evaluate` (同DiffDock框架) |
| 输出格式 | 同DiffDock(PDB) |
| 4090 GPU部署需求 | 同DiffDock |

**特点：**
- 在DiffDock基础上引入离散潜变量增强连续扩散
- 由DiffDock原作者+MIT/NVIDIA团队开发
- MIT License

---

## 5. 2024-2025最新SBDD/分子生成论文和开源仓库

### 5.1 DrugFlow (最推荐的新模型)

已在[1.3]详述。ICLR 2025, DiffSBDD作者后续模型, 当前SOTA。

### 5.2 Apo2Mol (最实用的新模型)

已在[3.3]详述。AAAI 2026, 首个Apo口袋条件生成模型。

### 5.3 BindDM

已在[2.3]详述。AAAI 2024, 子复合体增强的SBDD扩散模型。

### 5.4 GeoLDM (通用3D分子生成, 非口袋条件)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/MinkaiXu/GeoLDM |
| Stars | 278 |
| 论文标题 | Geometric Latent Diffusion Models for 3D Molecule Generation |
| 发表年份/会议 | ICML 2023 (arXiv: 2305.01140) |
| 预训练权重 | https://drive.google.com/drive/folders/1EQ9koVx-GA98kaKBS8MZ_jJ8g4YhdKsL |
| 推理脚本 | `eval_sample.py`, `eval_analyze.py` |
| 输出格式 | 内部格式(3D坐标), 霓自行转SDF/SMILES |
| SDF→SMILES转换难度 | 高 (输出为原始坐标+原子类型, 需手动构建分子对象) |
| 4090 GPU部署需求 | 需较大GPU内存(EGNN全连接), 建议4090 24GB |

**备注：** GeoLDM是通用3D分子生成模型(QM9/GEOM-Drugs), **不是口袋条件生成**, 但架构可改编为SBDD。

### 5.5 NeuralPLexer (蛋白柔性预测, 辅助SBDD)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/zrqiao/NeuralPLexer |
| Stars | 329 |
| 论文标题 | NeuralPLexer: Flexible Protein Structure Prediction with Deep Learning |
| 备注 | 不是分子生成模型, 但可预测蛋白柔性构象, 辅助口袋定义 |

### 5.6 DynamicBind (动态对接)

| 项目 | 信息 |
|------|------|
| GitHub | https://github.com/luwei0917/DynamicBind |
| Stars | 297 |
| 论文标题 | DynamicBind: Predicting Ligand-Specific Protein-Ligand Dynamics |
| 备注 | 不是生成模型, 是动态对接/构象预测, 可作为SBDD后端评估 |

### 5.7 其他值得关注的2024-2025论文(暂无开源代码)

以下论文在arXiv发表但暂未找到开源仓库:

1. **PocketGen** (2024): 蛋白口袋生成模型
2. **LigandDiff** (2024): 基于扩散的配体生成
3. **FlowSBDD** (推测): Flow Matching版SBDD
4. **CrossDiff** (2024): 跨模态扩散分子生成
5. **FLAG-SBDD** (推测): 基于FLAG架构的SBDD模型

---

## 综合对比表

| 模型 | 年份/会议 | 推荐度 | 推理脚本 | 输出格式 | SDF→SMILES | 单卡4090 | 特色功能 |
|------|----------|--------|---------|---------|-----------|---------|---------|
| **DrugFlow** | ICLR 2025 | ★★★★★ | src/generate.py | SDF | 低 | 可 | Flow+Markov Bridge, OOD检测, 侧链柔性, PA对齐 |
| **DiffSBDD** | Nat Comp Sci 2024 | ★★★★ | generate_ligands.py | SDF | 低 | 可 | Inpainting, 分子优化, Colab |
| **Apo2Mol** | AAAI 2026 | ★★★★ | sample_split.py | 内部 | 中 | 可 | Apo口袋条件(实用!), HF数据集 |
| **TargetDiff** | ICLR 2023 | ★★★★ | sample_diffusion.py/sample_for_pocket.py | meta→SDF | 中 | 可 | 亲和力预测双功能 |
| **Pocket2Mol** | ICLR 2023 | ★★★ | sample.py/sample_for_pdb.py | 内部→SDF | 中 | 可 | 自回归高效采样, taskset加速 |
| **DecompDiff** | ICML 2023 | ★★★ | sample_diffusion_decomp.py | meta→SDF | 中 | 可 | 分解先验, Vina Dock SOTA |
| **BindDM** | AAAI 2024 | ★★ | sample.py | meta | 中 | 可 | 子复合体增强 |
| **Pocket2Mol-RL** | 2024 | ★★ | 内部 | 内部 | 中 | 可 | RL微调 |
| **DiffDock-L** | ICLR 2024 | ★★★★(对接) | inference | PDB | N/A | 可 | 对接SOTA, HF在线, GUI |
| **DiffDock-Pocket** | 2024? | ★★(对接) | inference.py | PDB | N/A | 可 | 侧链柔性, ESMFold适配 |
| **DisCo-DiffDock** | ICML 2024 | ★★(对接) | evaluate | PDB | N/A | 可 | 离散潜变量增强 |

---

## 部署建议(4090 GPU)

**最低部署方案(单卡4090 24GB):**
- DrugFlow + DiffDock-L 组合: DrugFlow生成分子 → DiffDock-L评分/对接
- 环境建议: 使用DrugFlow的Docker镜像 `igashov/drugflow:0.0.3` 或独立conda环境
- 注意DrugFlow和TargetDiff/DecompDiff的conda环境可能冲突(PyG版本不同), 建议分环境部署

**推荐生产流水线:**
1. DrugFlow(生成SDF) → RDKit(SDF→SMILES过滤) → DiffDock-L(对接评分) → Vina/Gnina(精细对接)

---

## SDF→SMILES转换通用方案

所有输出SDF的模型都可统一用以下Python代码转换:

```python
from rdkit import Chem
mol = Chem.MolFromMolFile('generated.sdf')
if mol is not None:
    smiles = Chem.MolToSmiles(mol)
    print(smiles)
```

对于输出meta/内部格式的模型(TargetDiff, DecompDiff, Pocket2Mol), 需先运行其evaluate脚本, 或从3D坐标+原子类型手动构建RDKit Mol对象, 转换难度中等。
