# 代码演进日志 — Round 2

## 假设ID: H012

**假设**: 添加联芳基 C-C 键 Suzuki 逆合成规则可降低 trivial route 比例。

## 修改文件

### `src/synthesis_v2.py`

**修改内容**: 在 RETRO_RULES 列表最前面添加联芳基 C-C 键裂解规则

```python
# Suzuki 逆反应: Ar-Ar' → Ar-Br + (HO)2B-Ar'
("[c;R][c;R]", "[c:1][c:2]>>[c:1]Br.[c:2]B(O)O"),
```

**为什么放在最前面**: 联芳基断键应优先匹配，避免被后续的通用规则（如胺、醚）误匹配。

**安全性**:
- 仅匹配两个都在环中的芳香碳 (c;R)，避免误匹配链状 C-C 键
- 对稠环中的 C-C 键，产物无效（SMILES 验证自然拒绝），规则会 fallback

## 文献依据

- **Coscientist** (Boiko et al., 2023): 成功展示 Agent 自主规划 Suzuki 反应
- **LARC** (Baker et al., 2025): 规则覆盖率是逆合成质量的关键
- **ChemCrow** (Bran et al., 2024): 工具/库的丰富度直接决定 Agent 能力边界

## 验证

小测试确认对 Round 1 的 trivial biaryl 分子：
- 输入: `c1ccc(-c2ccc3ccccc3c2)cc1`
- 输出: `Brc1ccccc1.OB(O)c1ccc2ccccc2c1>>c1ccc(-c2ccc3ccccc3c2)cc1`
- Trivial: False ✅

## 代码自检

- `python3 -m py_compile src/synthesis_v2.py` ✅
