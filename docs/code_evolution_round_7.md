# 代码演进记录 — Round 7

## H015: 饱和含氮杂环断键规则 + _is_simple_molecule 修复

### 假设ID
H015

### 修改文件
`src/synthesis_v2.py` — `_is_simple_molecule` 函数 + `RETRO_RULES` 列表

### 修改内容

#### 1. _is_simple_molecule 修复（关键改动）
- **问题**: 多环芳烃（如邻甲基联苯，13 原子含苯环）被误判为"简单起始原料"
- **根因**: 旧逻辑仅检查 `atoms <= 15 && 含简单骨架`，不考虑环数
- **修复**: 新增环数检查 — 单环 ≤15 原子 → 简单；多环 ≤10 原子 → 简单；否则非简单
- **效果**: 联苯类分子正确进入 Suzuki 断键

#### 2. 新增 6 条逆合成规则
- 四氢异喹啉 (THIQ) → 逆 Pictet-Spengler
- 吲哚啉 → 逆还原环化
- 四氢喹啉 → 逆还原环化
- 饱和环 C-N 键 → 逆还原胺化
- 苄位 C-N 键 → 逆 Pictet-Spengler 变体
- 苯并氮杂环 → 逆 Bischler-Napieralski

### 文献依据
- LARC (Baker et al., 2025): 规则覆盖率 + 起始原料合理性
- Deep Lead Optimization (JACS 2024): 稠环体系的断键策略

### 验证状态
✅ 验证成功 — trivial ratio 0/10
