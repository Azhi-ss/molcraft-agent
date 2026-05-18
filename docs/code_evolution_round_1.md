# 代码演进日志 — Round 1 (H012)

## 假设 H012: 逆合成路线质量评分 + 规则优化

### 修改文件

#### 1. `src/synthesis_v2.py`

**新增函数 `score_route_quality()`**：
- 评估逆合成路线的化学合理性
- 评分维度：多步路线奖励 (+0.15)、多反应物路线 (+0.2)、复杂度比惩罚 (ratio > 3: -0.3)、单反应物检查
- 文献依据：LARC (Baker et al., 2025) Agent-as-a-Judge 框架

**修复吡啶逆合成规则**：
- SMARTS: `c1ccncc1` → `[c;R1]1[c;R1][c;R1][n;R1][c;R1][c;R1]1`
- `R1` 限制仅匹配孤立吡啶环（每个原子仅在 1 个 SSSR 环中）
- 排除喹啉、异喹啉、萘啶等稠环体系
- 验证通过：喹啉/异喹啉/萘啶 → False，吡啶 → True

**新增 6 条逆合成规则**：
1. 喹啉 → Friedländer 合成逆反应
2. 异喹啉 → Pictet-Spengler 逆反应
3. 喹唑啉 → 邻氨基苯甲酰胺 + 甲酸
4. 1,2,3-三唑 → Click Chemistry 逆反应
5. 肼/联氨 → 重氮还原

#### 2. `tools/pipeline.py`

**最终选择逻辑改为复合评分**：
- 候选池扩大至 3×n_top（30 个分子）
- 对每个候选分子：规划合成 + 路线质量评分
- 复合评分：`0.8 × BE_norm + 0.2 × route_quality`
- 按复合评分重新排序后选择 top N

### 文献依据
- LARC (Baker et al., 2025): Agent-as-a-Judge 评审路线可行性
- ChemCrow (Bran et al., 2024): 规则质量决定工具输出可靠性
- MOOSE-Chem (Yang et al., 2025): 多目标优化 — 合成性应为适应度维度之一

### 编译验证
- `src/synthesis_v2.py`: ✅ 编译通过
- `tools/pipeline.py`: ✅ 编译通过
- 函数测试: ✅ score_route_quality 正确区分好/坏路线
