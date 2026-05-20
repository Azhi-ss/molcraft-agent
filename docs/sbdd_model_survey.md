# 靶点条件3D分子生成模型开源仓库综合调研报告

> 调研时间：2026-05-20
> 调研目标：寻找能替代RDKit骨架变异方案的靶点条件3D分子生成模型，能在4090本地部署，接受PDB蛋白口袋输入，输出可转SMILES的3D分子

---

## 一、核心模型评估总表

| 模型 | GitHub链接 | Stars | 最近Push | 预训练权重 | PDB输入 | 输出→SMILES | 4090部署 | 推荐度 |
|------|-----------|-------|---------|-----------|---------|------------|---------|--------|
| **DiffSBDD** | arneschneuing/DiffSBDD | 504 | 2025-06 | ✅ Zenodo 8个ckpt | ✅ 直接--pdbfile | ✅ 输出SDF→SMILES | ⭐⭐⭐⭐⭐ | **首选** |
| **DrugFlow** | LPDI-EPFL/DrugFlow | 158 | 2026-01 | ✅ Zenodo 4个ckpt | ✅ --protein PDB+SDF参考配体 | ✅ SDF→SMILES | ⭐⭐⭐⭐⭐ | **强烈推荐** |
| **TargetDiff** | guanjq/targetdiff | 339 | 2024-01 | ✅ Google Drive | ✅ sample_for_pocket.py | ✅ reconstruct→SDF→SMILES | ⭐⭐⭐⭐ | 推荐 |
| **Pocket2Mol** | pengxingang/Pocket2Mol | 401 | 2023-11 | ✅ Google Drive | ✅ sample_for_pdb.py+中心坐标 | ✅ 直接输出SMILES+SDF | ⭐⭐⭐⭐ | 推荐 |
| **DecompDiff** | bytedance/DecompDiff | 74 | 2023-07 | ✅ Google Drive | ⚠️ 需lmdb数据集,无独立pdb入口 | ✅ 类似TargetDiff流程 | ⭐⭐⭐ | 可选 |
| **GraphBP** | divelab/GraphBP | 112 | 2023-07 (已迁入AIRS) | ✅ 仓库内model_33.pth | ✅ 基于CrossDocked数据 | ⚠️ 需后处理eval脚本 | ⭐⭐ | 较旧 |
| **IRDiff** | YangLing0818/IRDiff | 38 | 2024-09 | ⚠️ 需手动下载pretrained_models | ✅ 依赖TargetDiff流程 | ✅ 类似TargetDiff | ⭐⭐⭐ | 可选 |
| **Apo2Mol** | AIDD-LiLab/Apo2Mol | 32 | 2026-05(新发布) | ✅ apo2mol_checkpoint.ckpt | ✅ Apo蛋白口袋(无需参考配体!) | ✅ (待验证) | ⭐⭐⭐⭐ | **值得关注** |
| **D3FG** | 代码在EDAPINENUT/CBGBench | 5 | N/A | ⚠️ 需CBGBench框架 | ✅ 功能团级别生成 | ✅ 可转SMILES | ⭐⭐ | 较复杂 |
| **DrugDiff** | MarieOestreich/DrugDiff | 52 | 2025-03 | ✅ Zenodo | ❌ 非靶点条件(SELFIES) | ✅ SELFIES→SMILES | ❌ | **排除** |
| **DiffLinker** | igashov/DiffLinker | 386 | 2024+ | ✅ Zenodo | ✅ 可条件化于蛋白口袋 | ✅ SDF输出 | ⭐⭐⭐ | 辅助工具 |

---

## 二、重点模型详细分析

### 1. DiffSBDD ⭐⭐⭐⭐⭐ (首选推荐)

- **GitHub**: https://github.com/arneschneuing/DiffSBDD
- **论文**: Nature Computational Science (2024), ICLR风格
- **Stars**: 504，社区活跃度高
- **最近commit**: 2025-06-25 (持续更新)
- **预训练权重**: Zenodo提供8个checkpoint (CrossDocked/MOAD, Ca/fullatom, cond/joint)
- **PyTorch依赖**: PyTorch 1.12+ / PyTorch Lightning / CUDA 10.2+（新环境yaml支持更新版本）
- **输入格式**: 
  - `--pdbfile <pdb>` 直接接受PDB蛋白文件
  - `--ref_ligand A:330` 或 `--ref_ligand <sdf>` 定义口袋位置
  - `--resi_list A:1 A:2 ...` 可指定残基列表(无需参考配体)
- **输出格式**: SDF文件 (RDKit可直接 `Chem.MolToSmiles`)
- **特色功能**: 
  - 支持从头设计(de novo)
  - 支持子结构修复(inpainting/scaffold elaboration) ← 对molcraft很有用
  - 支持分子优化(optimize.py, 进化算法) ← 可优化SA/QED
  - 有Google Colab notebook
  - 有Docker支持
- **4090部署评估**: 
  - 推理只需单GPU，模型大小适中
  - 采样100个分子约需1-3分钟
  - **完全可行**，是最易部署的模型

**molcraft集成方案**: DiffSBDD可直接作为molcraft-agent的分子生成模块：
```python
# 调用方式
python generate_ligands.py checkpoint.ckpt \
  --pdbfile pocket.pdb \
  --outfile output.sdf \
  --resi_list A:100 A:101 ... \
  --n_samples 20

# 输出SDF→SMILES转换
from rdkit import Chem
suppl = Chem.SDMolSupplier('output.sdf')
for mol in suppl:
    if mol: smiles = Chem.MolToSmiles(mol)
```

---

### 2. DrugFlow ⭐⭐⭐⭐⭐ (强烈推荐，ICLR 2025最新)

- **GitHub**: https://github.com/LPDI-EPFL/DrugFlow
- **论文**: ICLR 2025 (Multi-domain Distribution Learning)
- **Stars**: 158
- **最近commit**: 2026-01-26
- **预训练权重**: Zenodo提供 DrugFlow/FlexFlow/DrugFlow+OOD/DrugFlow-PA 4个ckpt
- **PyTorch依赖**: PyTorch 2.2.1 + CUDA 12.1 (现代环境!)
- **输入格式**: 
  - `--protein examples/kras.pdb` (PDB文件)
  - `--ref_ligand examples/kras_ref_ligand.sdf` (SDF参考配体定义口袋)
  - `--pocket_distance_cutoff 8.0` (口袋半径)
- **输出格式**: SDF文件，自动计算SMILES和metrics
- **特色功能**:
  - **FlexFlow**: 联合采样蛋白侧链角度 + 分子 (考虑口袋柔性!)
  - Preference Alignment (DPO-like优化)
  - OOD检测 (uncertainty estimate)
  - Docker容器可用: `docker pull igashov/drugflow:0.0.3`
  - 当前DiffSBDD作者推荐的新模型("You can also try DrugFlow")
- **4090部署评估**: 
  - PyTorch 2.2+CUDA12.1完美匹配4090
  - Docker部署极其方便
  - **完全可行**

**注意事项**: DrugFlow需要ref_ligand来定义口袋位置，不能仅靠PDB+残基列表。molcraft需先用AutoDock Vina或AlphaFold获取参考配体位置。

---

### 3. TargetDiff ⭐⭐⭐⭐ (推荐，ICLR 2023经典)

- **GitHub**: https://github.com/guanjq/targetdiff
- **论文**: ICLR 2023 (3D Equivariant Diffusion)
- **Stars**: 339
- **最近commit**: 2024-01-10 (不再活跃更新)
- **预训练权重**: Google Drive提供
- **PyTorch依赖**: PyTorch 1.13.1 + CUDA 11.6 + PyG 2.2.0
- **输入格式**: 
  - `--pdb_path pocket.pdb` (有专用sample_for_pocket.py!)
  - 自动从PDB提取口袋原子
- **输出格式**: 
  - 原始输出为.pt文件(原子坐标+类型)
  - 通过reconstruct模块转为RDKit Mol → SDF → SMILES
  - **注意**: reconstruct使用OpenBabel，偶有重建失败
- **特色功能**:
  - 可同时做亲和力预测(affinity prediction)
  - 有Vina Docking评估内建
- **4090部署评估**: 
  - CUDA 11.6在4090上需要兼容性处理(4090需要CUDA 12+)
  - 推理可行，但环境配置可能需要调整PyTorch版本
  - **可行但需注意CUDA兼容**

---

### 4. Pocket2Mol ⭐⭐⭐⭐ (推荐，ML 2022经典)

- **GitHub**: https://github.com/pengxingang/Pocket2Mol
- **论文**: ICML 2022
- **Stars**: 401
- **最近commit**: 2023-11-16 (不再活跃)
- **预训练权重**: Google Drive提供 pretrained.pt
- **PyTorch依赖**: PyTorch 1.10.1 + CUDA 11.3 + PyG 2.0.0
- **输入格式**: 
  - `--pdb_path ./example/4yhj.pdb` (PDB文件)
  - `--center "32.0,28.0,36.0"` (需要手动指定口袋中心坐标!)
  - `--bbox_size 23.0` (口袋大小)
- **输出格式**: 
  - **直接输出SMILES** + SDF文件 ← 最方便转SMILES
  - 内建SMILES.txt输出
- **特色功能**:
  - 自回归式逐步生成(非扩散),采样更快
  - 直接输出SMILES，无需额外重建步骤
  - `taskset -c 0` 单CPU采样更快
- **4090部署评估**: 
  - CUDA 11.3 → 4090兼容性问题
  - 推理很轻量，内存占用小
  - **可行但需环境调整**

**关键限制**: 需要手动指定口袋中心坐标(3D坐标数值)，这在molcraft自动化流程中需要额外步骤来确定中心。

---

### 5. Apo2Mol ⭐⭐⭐⭐ (值得关注，AAAI 2026新模型)

- **GitHub**: https://github.com/AIDD-LiLab/Apo2Mol
- **论文**: AAAI 2026 (2026年5月刚发布代码)
- **Stars**: 32
- **预训练权重**: ✅ apo2mol_checkpoint.ckpt
- **核心创新**: **基于Apo蛋白口袋(无需参考配体/holo结构!)** ← 这对molcraft非常重要
  - 联合生成配体+口袋holo构象
  - 更符合实际应用场景(靶点蛋白通常没有已知配体)
- **输入格式**: Apo蛋白PDB文件(无需参考配体)
- **4090部署评估**: 待验证，但使用PyTorch Lightning+Hydra现代框架

---

### 6. DecompDiff ⭐⭐⭐ (可选，字节跳动)

- **GitHub**: https://github.com/bytedance/DecompDiff
- **论文**: NeurIPS 2023 workshop?
- **Stars**: 74
- **最近commit**: 2023-07 (不活跃)
- **预训练权重**: Google Drive提供
- **特色**: 分解先验(decomposed priors)，Vina Dock分数优于TargetDiff
- **限制**: 无独立sample_for_pdb脚本，需lmdb数据集格式
- **4090部署**: 类似TargetDiff环境，可行但需适配

---

### 7. GraphBP ⭐⭐ (较旧)

- **GitHub**: https://github.com/divelab/GraphBP (已迁入divelab/AIRS)
- **论文**: ICML 2022
- **Stars**: 112
- **最近commit**: 2023-07 (已停止维护)
- **预训练权重**: 仓库内model_33.pth (492KB，很小)
- **特色**: 自回归流模型，非扩散
- **限制**: 需48GB GPU训练(推理可行)，环境老旧(PyG 1.7.2)
- **4090部署**: 推理可行但环境老旧，不推荐优先

---

### 8. IRDiff ⭐⭐⭐ (可选)

- **GitHub**: https://github.com/YangLing0818/IRDiff
- **论文**: ICML 2024
- **特色**: 检索增强(retrieval-augmented)扩散模型
- **限制**: 依赖TargetDiff框架+PMINet预训练，安装复杂

---

### 9. DiffLinker ⭐⭐⭐ (辅助工具)

- **GitHub**: https://github.com/igashov/DiffLinker
- **论文**: Nature Machine Intelligence 2024
- **Stars**: 386
- **特色**: 片段连接(fragment linking) + 可条件化于蛋白口袋
- **用途**: 不是从头生成，而是连接已有片段。可作为molcraft的辅助模块。

---

### 10. DrugDiff ❌ (排除)

- **GitHub**: https://github.com/MarieOestreich/DrugDiff
- **原因**: 非靶点条件模型，基于SELFIES的潜在扩散，不接受蛋白口袋输入
- **用途**: 属性导向分子生成(QED/SA)，但不基于蛋白结构

---

## 三、Bohrium.com平台情况

- **平台地址**: https://www.bohrium.com (深势科技/DP Technology的AI for Science平台)
- **平台状态**: 正常运行，提供学术搜索、Notebook、Apps等功能
- **分子生成相关**: 
  - 有chem-sn化学版块(返回500错误，可能需登录)
  - Apps目录中发现: dfrocket(动力学), matmodeler(材料), matdesign(材料设计)
  - **未发现**专门的3D分子生成/SBDD应用
  - 平台主要侧重物理/材料模拟(DeepMD风格)
- **结论**: bohrium.com目前**没有**靶点条件3D分子生成模型的可直接应用，但可利用其GPU算力运行自己部署的模型

---

## 四、4090部署关键考量

### CUDA兼容性
- 4090需要 CUDA 12.x
- DrugFlow: PyTorch 2.2+CUDA12.1 → **完美匹配**
- DiffSBDD: 可适配CUDA12+PyTorch Lightning
- TargetDiff/Pocket2Mol: 需升级PyTorch到2.x并适配PyG

### 推理资源需求
| 模型 | GPU显存(推理) | 单次采样时间(100分子) | 难度 |
|------|------------|-------------------|------|
| DiffSBDD | ~2-4GB | ~1-3min | 低 |
| DrugFlow | ~4-6GB | ~2-5min | 低(Docker) |
| TargetDiff | ~3-5GB | ~5-10min | 中 |
| Pocket2Mol | ~2GB | ~1-2min | 中(需调环境) |

### molcraft-agent集成路径推荐

**方案A: DiffSBDD + DiffLinker组合 (最推荐)**
1. DiffSBDD做从头生成(de novo): 输入PDB口袋→输出SDF→SMILES
2. DiffSBDD做骨架修复(inpainting): 输入PDB+固定子结构→输出优化分子
3. DiffLinker做片段连接: 如果molcraft需要连接两个片段

**方案B: DrugFlow (最新最强)**
1. DrugFlow做从头生成: 输入PDB+参考配体SDF→输出SDF
2. FlexFlow可考虑口袋柔性
3. Docker部署最方便

**方案C: TargetDiff (经典稳定)**
1. 专用sample_for_pocket.py直接接受PDB
2. 需适配CUDA12环境
3. reconstruct模块偶有失败率

---

## 五、SMILES转换可行性

所有靶点条件3D分子生成模型的输出→SMILES路径:

1. **DiffSBDD/DrugFlow**: 直接输出SDF → RDKit `Chem.SDMolSupplier` + `Chem.MolToSmiles` → ✅简单
2. **TargetDiff**: reconstruct模块(OpenBabel+RDKit)重建分子 → 偶有失败 → ⚠️需要异常处理
3. **Pocket2Mol**: 内建直接输出SMILES.txt → ✅最简单
4. **通用**: 所有模型输出原子坐标+类型 → RDKit/OpenBabel重建 → 可行但有失败率

**关键问题**: 3D原子坐标→分子图重建(reconstruct)是所有扩散模型的共同难点。
- 扩散模型输出的是原子位置和类型，不保证化学合理性
- 需要从3D坐标推断键连接(基于距离阈值)
- 重建成功率通常在70-90%
- molcraft需要处理重建失败的case(返回None分子)

---

## 六、最终推荐

### 首选部署顺序:
1. **DiffSBDD** → 最成熟、最易用、PDB直接输入、SDF直接输出、有inpainting功能
2. **DrugFlow** → 最新ICLR 2025、PyTorch 2.2完美匹配4090、Docker部署、FlexFlow柔性口袋
3. **TargetDiff** → 经典ICLR 2023、有sample_for_pocket.py、但需环境适配
4. **Apo2Mol** → AAAI 2026、不需参考配体(Apo口袋)、但刚发布待验证

### 对molcraft-agent的具体建议:
- 当前binding_score是短板 → 3D生成模型直接在蛋白口袋中生成分子，天然具有更高亲和力
- 建议先部署DiffSBDD作为生成模块，替换RDKit骨架变异
- DiffSBDD的inpainting功能可用于保留关键子结构的同时优化其余部分
- 生成后仍需用AutoDock Vina验证亲和力(或使用TargetDiff的亲和力预测功能)

