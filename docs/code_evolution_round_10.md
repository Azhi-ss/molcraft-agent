# 代码演进报告 — Round 10

> 假设ID: H018 — BRICS 片段重组策略

---

## 修改摘要

| 项目 | 详情 |
|------|------|
| 假设ID | H018 |
| 修改文件 | `src/generator.py` |
| 文献依据 | MOOSE-Chem (2025), MolLEO (2024b), JACS 2024 |
| 改动行数 | ~180 行新增 + ~80 行修改 |

---

## 具体修改

### 1. 新增函数 `_brics_decompose_pool()` (generator.py:639-664)
- **功能**: 对分子池做 BRICS 分解，收集唯一片段
- **输入**: SMILES 列表
- **输出**: 含 attachment points 的 BRICS 片段 SMILES 列表
- **过滤**: 排除 < 2 重原子的片段

### 2. 新增函数 `_brics_recombine()` (generator.py:667-720)
- **功能**: 从 BRICS 片段库随机重组生成新分子
- **核心**: `BRICS.BRICSBuild(selected_frags)` → 取有效分子 → 类药过滤
- **参数**: 每次随机 2-4 个片段，每轮取 ≤3 个有效分子

### 3. 新增函数 `_build_fragment_db_from_scaffolds()` (generator.py:723-748)
- **功能**: 从 SCAFFOLDS 库的变异产物中提取 BRICS 片段（回退方案）
- **输出**: 97 个 BRICS 片段

### 4. 修复 `generate_molecules()` 的 combine 策略 (generator.py:751-855)
- **之前**: `"".join(parts)` 字符串拼接 SMILES
- **之后**: `_brics_decompose_pool()` → `_brics_recombine()` BRICS 片段重组
- **新增参数**: `fragment_pool` — 可选的外部分子池用于提取片段
- **回退**: combine 产量不足时自动补充 mutate 分子

### 5. 增强 `generate_with_docking_guidance()` (generator.py:944-1100)
- **gen=0**: combine 策略使用 BRICS 重组 + mutate 补充
- **gen>0**: 三种进化算子
  - 变异 (60%): 不变
  - Crossover (25%): H018 从 15% 提升到 25%
  - **BRICS 重组 (15%)**: 新增 — 从种子池 BRICS 分解→片段重组
- **文献依据**: MOOSE-Chem "recombination from population is essential"

---

## 验证状态

- [x] Python 编译通过
- [x] 导入测试通过 (`from generator import _brics_*` 成功)
- [x] BRICS 分解测试通过 (2 分子 → 5 片段)
- [x] BRICS 重组测试通过 (5 片段 → 6 有效分子)
- [x] 骨架库片段构建通过 (SCAFFOLDS → 97 片段)
- [ ] 完整 pipeline 运行（阶段四验证）
