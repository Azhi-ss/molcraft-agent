# 代码演进记录 Round 1

## 假设ID
H001: 收紧 SA score 过滤阈值

## 修改目标
将 `passes_filters()` 中的 SA score 阈值从 8.0 收紧至 6.0，并添加环数上限（MAX_RINGS=7）。
同时优化 `estimate_sa_score()` 的计算公式，使其更准确反映合成难度。

## 文献依据
- Deep Lead Optimization (JACS, 2024): 优质先导化合物的 SA score 通常在 2-5 之间
- LARC (Baker et al., 2025): 规则覆盖率和合成可行性是评估分子的核心维度
- 综述第3.2节: 分子生成需要约束在"可合成"空间内，避免浪费计算资源

## 修改文件
1. `src/evaluator.py` — 修改 `passes_filters()` 和 `estimate_sa_score()`

## 修改内容详情

### src/evaluator.py
**1. `estimate_sa_score()` 公式优化：**
- 环贡献因子: 0.5 → 1.0（多环结构显著增加合成难度）
- 添加稠环惩罚: 每个额外环系 +0.3（fused ring penalty）
- 螺环惩罚: 1.0 → 1.5（螺环形成挑战性高）
- 桥头原子惩罚: 1.5 → 2.0（桥头结构极难合成）
- 保持手性中心和可旋转键贡献不变

**2. `passes_filters()` 新增参数：**
- `max_sa=6.0`: 从 8.0 降低至 6.0
- `max_rings=7`: 限制分子中环数不超过 7 个

## 代码自检
```python
python3 -m py_compile src/evaluator.py
```
