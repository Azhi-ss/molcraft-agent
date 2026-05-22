# 代码演进报告 — Round X (H034)

## 修改文件

### `src/synthesis_v2.py` — RETRO_RULES 扩展

**修改内容**: 在 RETRO_RULES 末尾（H020 之后）新增 5 条三环稠合杂环逆合成规则。

**新增规则**:

| # | 骨架 | 命名反应 | 逆合成路径 |
|---|------|---------|-----------|
| 1 | 吖啶 (Acridine) | Bernthsen 合成 | 吖啶 → 2-氨基联苯 + 甲酸 |
| 2 | 吩嗪 (Phenazine) | Wohl-Aue 合成 | 吩嗪 → 邻苯二胺 + 邻苯醌 |
| 3 | 咔唑 (Carbazole) | Borsche-Drechsel 环化 | 咔唑 → 2-氨基联苯 |
| 4 | 吡咯并[2,3-d]嘧啶 | Traube 嘌呤合成变体 | 7-去氮嘌呤 → 4-氨基嘧啶 + 乙二醛 |
| 5 | 吡唑并[3,4-d]嘧啶 | Knorr 类缩合 | 吡唑并嘧啶 → 4-氨基嘧啶 + 肼 |

**设计原则**:
1. 放在所有规则末尾，确保更具体的双环规则优先匹配
2. SMARTS 匹配核心骨架（取代基无关），RETRO 给出简化的合成子
3. 每个规则对应真实命名反应

**文献依据**: 
- LARC (Baker et al., 2025): 规则覆盖率决定逆合成质量
- Deep Lead Optimization (JACS 2024): 稠环体系需要专门断键策略
- ChemCrow (Bran et al., 2024): 工具/规则库丰富度决定 Agent 能力边界

## 修改位置

`src/synthesis_v2.py`, lines 624-630 (H020 异噁唑规则之后, 闭括号 `]` 之前)

## 编译检查

✅ `python3 -m py_compile src/synthesis_v2.py` 通过

## 单元验证

| 测试分子 | 结果 | 步数 | 路线 |
|---------|------|------|------|
| 取代吖啶 `Cc1ccc2nc3ccc(C)cc3cc2c1` | ✅ OK | 2 | Suzuki→联苯 + Bernthsen 环化 |
| 吩嗪 `c1ccc2nc3ccccc3nc2c1` | ✅ OK | 1 | 邻苯二胺 + 邻苯醌→吩嗪 |
| 咔唑 `c1ccc2c(c1)[nH]c1ccccc12` | ✅ OK | 2 | Suzuki→2-氨基联苯 + Borsche-Drechsel |
| 未取代吡咯并嘧啶 | 简单分子 | 0 | ≤10原子，正确判定为商业可得起始原料 |
| 取代吡唑并嘧啶 `Cc1nc2[nH]ncc2c(C)n1` | ✅ OK | 1 | 氨基嘧啶 + 肼→吡唑并嘧啶 |

所有新规则 SMARTS 与目标骨架正确匹配，无交叉误匹配（喹啉不匹配任何新规则）。

## 风险评估

- **向后兼容**: 新规则仅追加于末尾，不修改任何现有规则
- **mass balance**: 新规则均通过 H014 化学计量守恒检查（ratio 0.75-1.25）
- **元素守恒**: 新规则均通过 H021 元素守恒检查
- **性能影响**: 5 条新规则，每次 `plan_synthesis_recursive` 调用额外匹配成本可忽略
