# 代码演进报告 — Round 22 (H026)

## 假设 H026: 扩展大芳香骨架库

### 改动文件
- `src/generator.py` — SCAFFOLDS 列表

### 改动内容
新增 14 个大芳香/杂芳烃骨架，同时修复 5 个已有无效 SMILES：

**新增骨架 (H026):**
| 骨架 | 类型 | 环数 | 重原子 |
|------|------|------|--------|
| 蒽 (anthracene) | 3环线性芳烃 | 3 | 14 |
| 菲 (phenanthrene) | 3环角形芳烃 | 3 | 14 |
| 芘 (pyrene) | 4环芳烃 | 4 | 16 |
| 芴 (fluorene) | 联苯+桥 | 3 | 13 |
| 荧蒽 (fluoranthene) | 4环芳烃 | 4 | 16 |
| 苯并[a]蒽 | 3环芳烃 | 3 | 14 |
| 吖啶 (acridine) | N代蒽 | 3 | 14 |
| 吩嗪 (phenazine) | 二氮代蒽 | 3 | 14 |
| 吩噻嗪 (phenothiazine) | S+N三环 | 3 | 14 |
| 咔唑 (carbazole) | N代芴 | 3 | 13 |
| 二苯并呋喃 (dibenzofuran) | O三环 | 3 | 13 |
| 二苯并噻吩 (dibenzothiophene) | S三环 | 3 | 13 |
| 苯并[h]喹啉 | N杂三环 | 3 | 14 |
| 喹喔啉 (quinoxaline) | 苯并吡嗪 | 2 | 10 |

**修复已有无效 SMILES:**
- 喹唑啉: `c1ccc2c(c1)Ncnc2` → `c1ccc2cncnc2c1`
- 苯并咪唑: `c1ccc2c(c1)ncn2` → `c1ccc2nc[nH]c2c1`
- 吲哚啉: `c1ccc2c(c1)[nH]c2` → `c1ccc2CCNc2c1`
- 嘌呤: `c1ncnc2c1ncn2` → `c1[nH]cnc2ncnc12`
- 吲哚嗪: `c1ccc2nccc2c1` → `c1ccc2cccn2c1`

### 文献依据
- Deep Lead Optimization (JACS 2024) §3.2: 骨架多样性决定化学空间探索效率；"the variety of parent skeletons should be notably broad"
- H025 (REJECTED): Vina 惩罚极性取代基 → 偏好疏水体系
- H027 (VERIFIED): 放宽过滤器后 BE +4.2% → 更大分子有优势

### 代码自检
- ✅ 所有 84 个骨架 SMILES 可解析
- ✅ generator.py 编译通过
- ✅ 全模块 import 链正常

### 改动原则
- 仅新增骨架 SMILES，不改动任何逻辑代码
- 保留现有接口不变
- 不影响其他模块
