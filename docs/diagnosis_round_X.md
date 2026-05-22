# 诊断报告 — Round X (2026-05-22 Session)

## 1. 禁重检查

已检查 knowledge_base.md 中所有 REJECTED 条目:
- ❌ H019: Suzuki Byproduct Atom Balance — trivial ratio regressed (不可重提)
- ❌ H025: Expanded Substituent Library (Polar Pharmacophores) — BE degraded (不可重提)
- ❌ H026: Expanded Large Aromatic Scaffold Library — BE degraded (不可重提)

新假设方向与上述 REJECTED 条目均不重叠。

## 2. 外部新知识搜索

### LKM 搜索
搜索 TYK2 inhibitor + docking 相关文献，获得以下关键发现：
- TYK2 选择性抑制剂的系统性虚拟筛选仍是不足领域（gcn_3b95e74852f24356）
- AutoDock Vina 评分函数的 docking score 差异 ~1 kcal/mol 对应显著亲和力变化（gcn_83f59ea641254599）
- DecompOpt (2024): 可控分解扩散模型用于基于结构的分子优化（arXiv:2403.13829）

### arXiv 搜索
- RetroReasoner (2026-03): 基于推理 LLM 的逆合成预测，使用 RL 优化（arXiv:2603.12666）
- Margin-calibrated Classifier Guidance (2026-05): 属性驱动的合成规划（arXiv:2605.13101）

**可操作性评估**: 最新论文以深度学习为主，与项目规则基础架构差距大。直接可落地策略有限。

## 3. 源码审查

已审查文件：
- `src/synthesis_v2.py`: RETRO_RULES 已覆盖 50+ 条规则（Suzuki, amide, sulfonamide, ester, ether, amine, heterocycles, Diels-Alder, lactone, epoxide, pyrazole, sp3 C-C 等）
- `tools/pipeline.py`: H032 后过滤器已实现（line 321-336），排除 smiles>>smiles trivial 路线
- `src/generator.py`: 使用 SCAFFOLDS + KINASE_HINGE_SCAFFOLDS，H021 激酶偏置已启用

H032 代码逻辑正确：
```python
valid_candidates = [c for c in scored_candidates if c["route"] != f"{c['mol_smiles']}>>{c['mol_smiles']}"]
```
上一会话因工具模块缓存导致未实测（session progress 记载）。

## 4. 瓶颈诊断

### 当前基线 (H031 VERIFIED)
| Metric | Value |
|--------|-------|
| Best BE | -9.902 |
| Avg BE | -8.642 |
| Trivial ratio | 2/10 |

### 瓶颈识别

**瓶颈 1: H032 未实测验证**
- 代码正确但未被 live pipeline 测试
- 预期效果: 排除 smiles>>smiles 分子，trivial ratio 从 2/10 → ≤1/10

**瓶颈 2: 四环稠合体系无断键规则**
- 知识库记载: "2 trivial are smiles>>smiles type (tetracyclic scaffolds w/o synthesis rules)"
- 当前 RETRO_RULES 覆盖了常见双环/三环杂环（喹啉、异喹啉、吲哚、苯并呋喃等）
- 但三环/四环稠合杂环体系（如 acridine, phenazine, carbazole, pyrrolopyrimidine）仍无匹配规则
- 这些体系在激酶抑制剂中高频出现（hinge-binding scaffolds）

**瓶颈 3: BE 提升空间有限**
- 当前最佳 BE -9.902 已处于 Vina 评分函数天花板附近
- Vina 偏向疏水/扁平芳环体系（H025 结论），进一步优化需谨慎

### 工具天花板评估
- Vina 对接：已有 consensus docking (H009)、docking guidance (H002)
- 逆合成：规则覆盖率已达 50+，但三/四环稠合体系为明确缺口
- 分子生成：SCAFFOLDS 库丰富（H022），变异策略多样

## 5. 假设提出

### H034: Tricyclic Fused Heterocycle Retrosynthesis Rules + H032 Live Verification

**来源**: LARC (Baker et al., 2025) — 规则覆盖率决定逆合成质量

**问题**: 当前 RETRO_RULES 缺少以下在激酶抑制剂中高频出现的稠合杂环断键规则：
1. **Acridine** (二苯并[b,e]吡啶): 三环含氮芳环 → 逆 Ullmann/环化
2. **Phenazine** (二苯并[b,e]吡嗪): 三环双氮芳环 → 逆缩合
3. **Carbazole** (二苯并[b,f]吡咯): 三环含氮芳环 → 逆 Cadogan/Borsche-Drechsel
4. **Pyrrolo[2,3-d]pyrimidine** (7-deazapurine): 激酶 hinge-binder 核心
5. **Pyrazolo[3,4-d]pyrimidine**: 常见激酶抑制剂骨架

**方案**:
1. 在 RETRO_RULES 中新增 5 条逆合成规则，覆盖上述稠合杂环
2. 每条规则对应真实命名反应或仿生合成路径
3. 同时运行 pipeline 实测 H032 + H034 效果

**预期效果**: 
- Trivial ratio: 2/10 → 0/10
- Best BE: -9.902 ± 0.3 (不应显著退化)
- 新增规则覆盖的分子应获得 ≥1 step 非平凡路线

**验证标准**:
- H032 后过滤器正确排除 smiles>>smiles 路线
- 新规则至少为 1 个分子提供非平凡逆合成路线
- Best BE 不低于 -9.5 kcal/mol（保守门槛）
