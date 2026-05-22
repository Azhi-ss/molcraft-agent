# MolCraft Agent 🧬

> **靶向分子研发与合成规划智能体**
>
> 第四届世界科学智能大赛 —— AI4S 智能体 CNS 挑战赛（任务2）

## 项目简介

MolCraft Agent 是一个 **Agent 驱动** 的靶向药物小分子设计与合成路线规划系统。该系统利用大语言模型作为"大脑"，通过自定义化学工具自主完成从 **分子生成**、**分子对接** 到 **逆合成路线规划** 的全流程闭环，并支持基于对接结果的 **迭代优化**。

运行模式为 **Agent 模式**：LLM 自主决策，按「文献解析 → 瓶颈诊断 → 代码演进 → 实验验证」四阶段循环迭代。`result.log` 自动记录 Agent 完整决策过程，用于证明方案由智能体流程产生。

## 赛题背景

在新药研发中，先导化合物的发现与优化是耗时最久、投入最高的环节之一。本赛题要求构建高度自动化的智能体，在规定时间内生成符合要求的药物小分子，同时给出对应的合成路线。

**评分维度：**
- 与靶点的结合能（越低越好）
- 分子结构合理性
- 可合成性
- 起始原料的可及性
- 合成路线的经济性

## 快速开始

### 1. 环境准备

本项目使用 [uv](https://github.com/astral-sh/uv) 管理 Python 环境。

```bash
# 安装 uv（若尚未安装）
curl -LsSf https://astral.sh/uv/install.sh | sh

# 进入项目目录
cd molcraft-agent

# 同步创建虚拟环境并安装所有依赖（Python ≥ 3.12）
uv sync --python 3.12
source .venv/bin/activate
```

### 2. 运行方式

#### 方案 A：Pipeline 模式（直接生成结果）

如果你需要**直接产出** `result.csv` 和 `result.log`，不需要迭代优化，用这个模式：

```bash
source .venv/bin/activate

# 生成 50 个分子，取对接得分前 10 的做逆合成规划
python tools/pipeline.py --n-generate 50 --n-top 10
```

Pipeline 是纯脚本执行，不调用 LLM，走确定性流程：分子生成 → 对接筛选 → 逆合成路线规划。结果直接输出到 `output/`。

#### 方案 B：Agent 模式（迭代优化）

如果你希望 Agent **自主迭代**优化分子，提升对接得分和合成可行性，用这个模式：

```bash
source .venv/bin/activate

# 默认运行 1 次迭代（约 30-60 分钟），最长 90 分钟，每轮最多 1000 步
python main.py

# 迭代 3 次
python main.py --iterations 3

# 最长运行 60 分钟
python main.py --max-minutes 60

# 每轮最多 500 步（Kimi CLI 默认步数上限）
python main.py --max-steps 500
```

Agent 将读取 `program.md` 中的指令，自主执行：
1. **文献解析** — 读取 `papers/` 中的论文，提取关键架构和方法
2. **瓶颈诊断** — 分析现有代码，对比论文方法找出差距
3. **代码演进** — 根据最高优先级假设，自主修改代码
4. **实验验证** — 调用化学工具运行实验，量化评估改进效果
5. **科学迭代** — 验证成功则深化，失败则分析原因并转向新假设

Agent 产出将保存在 `docs/` 和 `output/` 目录中。

**运行控制**：通过 `--iterations`（默认 1）、`--max-minutes`（默认 90）和 `--max-steps`（默认 1000）三重限制，防止 Agent 无限运行。

**人类迭代策略**：编辑 `program.md` 即可调整 Agent 行为，无需修改代码。

### 3. 打包提交

```bash
python3 -c "import zipfile; z=zipfile.ZipFile('result.zip','w'); z.write('output/result.csv','result.csv'); z.write('output/result.log','result.log'); z.close()"
```

## 项目结构

```
molcraft-agent/
├── data/                       # 输入数据
│   ├── target.pdb              # 蛋白质靶点文件
│   └── receptor.pdbqt          # 准备好的受体文件（自动生成）
├── src/                        # 核心化学计算引擎
│   ├── config.py               # 配置模块（路径、对接参数、过滤条件）
│   ├── receptor.py             # 受体准备（PDB → PDBQT）
│   ├── generator.py            # 分子生成（RDKit 变异 / 对接引导生成）
│   ├── evaluator.py            # 分子性质评估（QED / Lipinski / SA score）
│   ├── docking.py              # 分子对接（AutoDock Vina）
│   └── synthesis_v2.py         # 逆合成（SMARTS 模板 + BRICS fallback）
├── molcraft_agent/             # Agent SDK 集成层
│   ├── __init__.py
│   ├── tools.py                # 4 个自定义化学工具定义
│   └── experiments.py          # 实验记录模块
├── tools/                      # CLI 工具
│   ├── pipeline.py             # 纯脚本入口（无 LLM，仅调试）
│   ├── generate.py             # 分子生成 CLI
│   ├── dock.py                 # 对接 CLI
│   ├── evaluate.py             # 评估 CLI
│   ├── synthesize.py           # 合成规划 CLI
│   └── git_advance.py          # Git 实验管理（advance / revert）
├── papers/                     # 参考论文（Agent 会读取并分析）
│   ├── autonomous_agents_survey.md
│   └── coscientist.md
├── docs/                       # Agent 产出文档
│   ├── literature_analysis_round_*.md
│   ├── diagnosis_round_*.md
│   ├── code_evolution_round_*.md
│   └── experiment_round_*.md
├── output/                     # 结果输出
├── tutorial/                   # 使用教程（新手从这里开始）
│   ├── 00-baseline.md          # 跑通第一个 Agent
│   ├── 01-understanding.md     # 理解 Agent 工作原理
│   └── 02-improvement.md       # 提升 Agent 效果
├── program.md                  # Agent 指令书（人机接口，人类编辑）
├── agent.yaml                  # Agent 工具配置
├── main.py                     # Agent 模式启动器
├── pyproject.toml              # 项目配置
├── uv.lock                     # 依赖锁定文件
├── .env.example                # 环境变量模板
└── README.md                   # 本文件
```

## 架构设计

本项目采用 **Program-driven** 的自主科研架构：

- **`program.md`** —— 人类编写的 Agent 指令书，Agent 严格遵循。这是**唯一的人机接口**，人类通过编辑它来迭代研究策略。
- **`experiments.jsonl`** —— 结构化实验记录，每次工具调用自动追加，完整可追溯。
- **`main.py`** —— 极简启动器，只负责读取 `program.md` 并启动 Agent，自身不携带业务逻辑。运行参数：`--iterations`（迭代次数，默认 1）、`--max-minutes`（最大运行时间，默认 90）。
- **四阶段闭环** —— Agent 自主循环：文献解析 → 瓶颈诊断 → 代码演进 → 实验验证 → 迭代优化。

```
┌─────────────────────────────────────────────────────────────┐
│  Human                                                    │
│  └── 编辑 program.md（迭代策略）                            │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│  Agent (LLM)                                              │
│  ├── 读取 program.md                                       │
│  ├── 阶段一：文献解析（papers/）                            │
│  ├── 阶段二：瓶颈诊断（src/ 代码分析）                       │
│  ├── 阶段三：代码演进（StrReplaceFile 修改）                │
│  ├── 阶段四：实验验证（调用化学工具）                      │
│  └── LOOP：验证成功则深化，失败则转向新假设                  │
└──────────────────────┬──────────────────────────────────────┘
                       │
          ┌────────────┼────────────┐
          ▼            ▼            ▼
  ┌──────────────┐ ┌──────────┐ ┌─────────────┐
  │ RDKit 变异   │ │ Vina 对接 │ │ SMARTS 模板 │
  │ 进化迭代     │ │ 结合能   │ │ 逆合成      │
  └──────────────┘ └──────────┘ └─────────────┘
```

### 分子生成
- 基于 27 种药物样骨架的随机变异（添加取代基、替换原子、插入连接子）
- 支持进化迭代：从对接成功的种子分子变异产生后代
- 支持对接引导生成：利用实时对接反馈指导生成方向
- QED、Lipinski 五规则、分子量（150–500）、LogP（−0.5–5）、SA score 多维度过滤

### 分子对接
- 受体：Meeko 准备 PDBQT
- 配体：RDKit 3D 构象 + Meeko 转换为 PDBQT
- 对接引擎：AutoDock Vina（exhaustiveness=8）

### 逆合成规划
- 基于 SMARTS 的反应模板切断规则
- 支持：酰胺、磺酰胺、酯、芳基醚、还原胺化、Buchwald-Hartwig、脲、硝基还原等
- 回退策略：BRICS 碎片化 → 平凡路线

## 自定义化学工具

Agent 通过以下工具与化学计算引擎交互：

| 工具 | 功能 | 参数 |
|------|------|------|
| `generate_molecules` | 生成候选分子 | `strategy`, `n_molecules` |
| `dock_molecules` | 批量对接 | `smiles_list` |
| `evaluate_molecule` | 评估药物性质 | `smiles` |
| `plan_synthesis` | 规划逆合成路线 | `smiles` |

工具定义在 `molcraft_agent/tools.py`，通过 `agent.yaml` 注册到 Kimi Agent Runtime。

Agent 还拥有 13+ 个内置工具（文件读写、代码修改、Shell 执行、Web 搜索、任务管理等），支持自主科研全流程。

## 环境变量配置

项目通过 `.env` 文件加载所有 API 密钥和服务地址。复制模板文件并填入真实值：

```bash
cp .env.example .env
# 编辑 .env 填入密钥
```

### 【必填】大语言模型（Agent 的大脑）

任意选择一个兼容 `/v1/chat/completions` 的 LLM 服务，取消 `.env` 中对应方案的注释：

| 方案 | Base URL | Key 获取 | 推荐模型 | 费用参考 |
|------|----------|---------|---------|---------|
| **DeepSeek** ✅推荐 | `https://api.deepseek.com/v1` | [platform.deepseek.com](https://platform.deepseek.com) → API Keys | `deepseek-v4-pro` | ¥30 可跑完整一轮 |
| Kimi | `https://api.moonshot.cn/v1` | [platform.moonshot.cn](https://platform.moonshot.cn) → API Keys | `kimi-k2.6` | 新用户有免费额度 |
| OpenAI | `https://api.openai.com/v1` | [platform.openai.com](https://platform.openai.com) → API Keys | `gpt-4o` | 需绑信用卡 |
| 阶跃星辰 | `https://api.stepfun.com/v1` | [platform.stepfun.com](https://platform.stepfun.com) → API Keys | `step-2-16k` | 免费额度 |
| 智谱 GLM | `https://api.z.ai/api/paas/v4` | [open.bigmodel.cn](https://open.bigmodel.cn) → API Keys | `glm-4.6` | 免费额度 |
| 本地模型 | `http://localhost:11434/v1` | 无需 key（填 `ollama`） | `qwen2.5:14b` | 免费 |

对应的环境变量（三行一组）：
```
LLM_BASE_URL=https://api.deepseek.com/v1
LLM_API_KEY=sk-your-real-api-key
LLM_MODEL=deepseek-v4-pro
```

### 【可选但推荐】扩散模型（GPU 服务器）

用于 PocketXMol 口袋感知分子生成，需在 GPU 服务器（如 4090）上单独部署服务。部署方式见 `docs/diffusion_server.md`。

```
DIFFUSION_API_URL=http://your-gpu-server:8001
```

### 【可选】Bohrium LKM 科学知识图谱

Agent 通过 `SearchLKM` 工具检索科学文献和 claims。获取 key：

1. 访问 [bohrium.dp.tech](https://bohrium.dp.tech) → 注册/登录
2. 工作台 → 个人设置 → API 令牌 → 创建
3. 复制 access key 填入 `.env`

```
LKM_ACCESS_KEY=your-lkm-access-key
```

不配置时 SearchLKM 工具会返回错误提示，不影响核心功能（分子生成、对接、逆合成）。

### 【可选】逆合成 API（IBM RXN）

```
RXN_API_KEY=your-rxn-api-key
```

获取地址：[rxn.res.ibm.com](https://rxn.res.ibm.com)。不配置时使用内置 SMARTS 模板，效果足够竞赛提交。

### 完整 .env 示例

```bash
# === 必填 ===
LLM_BASE_URL=https://api.deepseek.com/v1
LLM_API_KEY=sk-your-real-api-key
LLM_MODEL=deepseek-v4-pro

# === 可选 ===
DIFFUSION_API_URL=http://10.0.0.1:8001    # GPU 服务器
LKM_ACCESS_KEY=your-lkm-access-key         # Bohrium 知识图谱
RXN_API_KEY=                               # IBM RXN（不用可以不填）
```

## 核心运行参数

以下参数位于 `src/config.py`，按比赛规格（TYK2 / PDB 5C01）调优：

### 分子对接（TYK2 靶点）
| 参数 | 当前值 | 说明 |
|------|--------|------|
| `DOCKING_CENTER` | `[21.86, -0.41, 29.93]` | 对接盒子中心（几何口袋检测活性位点） |
| `DOCKING_SIZE` | `[28.0, 28.0, 28.0]` | 盒子大小（Å³=21952，低于 Vina 推荐上限 27000） |
| `DOCKING_EXHAUSTIVENESS` | `32` | 搜索强度，越高越可能发现深结合模式 |

### 类药性质过滤
| 参数 | 值 | 说明 |
|------|-----|------|
| `MAX_MW` | 550 | 最大分子量 |
| `MIN_MW` | 150 | 最小分子量 |
| `MAX_LOGP` | 5.0 | 最大脂溶性 |
| `MIN_LOGP` | -0.5 | 最小脂溶性 |
| `MAX_TPSA` | 140 | 最大极性表面积 |
| `MIN_QED` | 0.3 | 最低类药性 |

### 分子生成
| 参数 | 值 | 说明 |
|------|-----|------|
| `N_GENERATED_MOLECULES` | 50 | 每次生成数量 |
| `N_TOP_MOLECULES` | 10 | 保留最优数量 |

### 迭代控制
| CLI 参数 | 默认值 | 说明 |
|----------|--------|------|
| `--iterations` | 1 | 循环迭代次数 |
| `--max-minutes` | 90 | 单次最长运行时间 |
| `--max-steps` | 1000 | 每轮最大决策步数 |

## 依赖列表

- Python ≥ 3.12
- RDKit（分子操作）
- AutoDock Vina（分子对接）
- Meeko（PDBQT 准备）
- Open Babel（格式转换）
- kimi-agent-sdk（Agent 运行时）
- kimi-cli（工具框架）

完整依赖见 `pyproject.toml`。运行 `uv sync --python 3.12` 自动安装全部依赖并生成 `uv.lock`。

## 许可协议

MIT License

## 致谢

第四届世界科学智能大赛 —— AI4S 智能体 CNS 挑战赛
