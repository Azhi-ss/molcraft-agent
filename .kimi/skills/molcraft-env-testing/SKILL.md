---
name: molcraft-env-testing
description: 外部工具独立环境测试规范。当需要测试 FPocket、Uni-Dock、DeepSite 等不在项目依赖中的工具时读取。
---

# 独立环境工具测试规范

如果判断需要测试某个外部工具（如 FPocket、Uni-Dock、DeepSite 等），必须在**独立的临时 conda 环境**中执行，不能污染项目 `.venv/`。

## 操作流程

```bash
# 步骤 1：创建独立临时环境（用完即毁）
conda create -n molcraft-tools-$(date +%Y%m%d) -c conda-forge fpocket -y

# 步骤 2：在独立环境中运行工具（用 conda run，不需要手动 activate）
conda run -n molcraft-tools-$(date +%Y%m%d) fpocket -f data/target.pdb

# 步骤 3：销毁环境（用完必删）
conda env remove -n molcraft-tools-$(date +%Y%m%d) -y
```

## 重要规则

- 每次测试都新建环境，测试完立刻删除
- 如需长期保留（替换了工具），另建一个 `molcraft-tools` 固定环境
- 不要往 `.venv/` 或 conda base 里装任何工具测试相关的包
- 接口对齐：新模块输出格式与原有模块兼容（同样的返回值结构）
- 默认保留原有工具作为可选项，验证确实更好再切换