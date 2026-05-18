# 结构化 JSON Lines 日志系统设计文档

**日期**: 2026-05-17
**版本**: 1.0
**状态**: 待审核

## 1. 背景与目标

### 1.1 问题
- 原系统的 `result.log` 是纯文本/Markdown 格式，导致平台 LLM 评分不稳定
- 关键指标（结合能、SA、trivial ratio）分散在日志各处
- 格式不一致，每轮生成的报告结构有变化

### 1.2 平台要求
- `result.log` 内容必须是 JSON 格式（文件名不带 .json 后缀）
- 提交时需要打包成 zip: `zip result.zip result.csv result.log`
- 纯 JSON 评分成功率 100%，纯文本评分报错率高

### 1.3 目标
- 将 `result.log` 从纯文本改为 **JSON Lines 格式**
- Agent 运行过程中实时流式写入事件
- 运行结束后自动打包 `result.zip`
- 终端输出保持人类可读性，不影响实时监控

## 2. 数据格式设计

### 2.1 JSON Lines 格式
每行一个独立的 JSON 对象，支持流式追加：

```json
{"type": "start", "timestamp": "2026-05-17T15:13:25", "round": 3, "max_minutes": 90}
{"type": "stage", "name": "启动检查", "status": "completed", "duration_seconds": 12.5}
{"type": "docking_progress", "current": 10, "total": 60, "success_rate": 1.0}
{"type": "molecule", "smiles": "Cc1cc(C)c2ccccc2c1O", "binding_energy": -8.098, "sa_score": 2.7, "rings": 2, "synthesis_steps": 0, "trivial_route": true, "synthesis_route": "..."}
{"type": "hypothesis_validation", "hypothesis_id": "H001", "success": true, "conclusion": "ACCEPTED", "changes": {"sa_threshold": {"from": 8.0, "to": 6.0}, "max_rings": 7}}
{"type": "metrics", "top10_avg_be": -7.761, "best_be": -8.098, "trivial_ratio": 0.1, "avg_synthesis_steps": 1.11, "docking_success_rate": 1.0, "non_trivial_count": 9}
{"type": "end", "timestamp": "2026-05-17T15:20:31", "duration_seconds": 426, "status": "success"}
```

### 2.2 事件类型定义

| 事件类型 | 触发时机 | 关键字段 |
|----------|----------|----------|
| `start` | Agent 启动时 | timestamp, round, max_minutes, max_iterations |
| `stage` | 每个阶段开始/结束 | name, status, duration_seconds |
| `docking_progress` | 每完成一批分子对接 | current, total, success_rate |
| `molecule` | 每个分子的逆合成分析完成 | smiles, binding_energy, sa_score, rings, synthesis_steps, trivial_route, synthesis_route |
| `hypothesis_validation` | 假设验证完成 | hypothesis_id, success, conclusion, changes |
| `metrics` | 所有实验完成后汇总 | top10_avg_be, best_be, trivial_ratio, avg_synthesis_steps, docking_success_rate, non_trivial_count |
| `end` | Agent 正常结束或异常退出 | timestamp, duration_seconds, status, error_message |

## 3. 系统架构

### 3.1 核心组件

```
main.py
├── StructuredLogger (替换 TeeLogger)
│   ├── log_event(type: str, **data)    # 写入 JSON Lines 到文件
│   ├── _print_human_readable(event)    # 终端打印格式化文本
│   └── flush() / close()
│
├── Agent 执行 (kimi_agent_sdk.prompt)
│   ├── 启动 → log_event("start", ...)
│   ├── 阶段切换 → log_event("stage", ...)
│   ├── 对接完成 → log_event("docking_progress", ...)
│   ├── 分子合成 → log_event("molecule", ...)
│   └── 验证完成 → log_event("hypothesis_validation", ...)
│
└── finally 块
    ├── 写入 metrics 汇总
    ├── 写入 end 事件
    ├── 调用 zip_results() → result.zip
    └── 日志关闭
```

### 3.2 文件修改清单

| 文件 | 修改内容 |
|------|----------|
| **main.py** | 1. 将 `TeeLogger` 类替换为 `StructuredLogger`<br>2. 添加 `zip_results()` 函数<br>3. 在 finally 块调用打包逻辑 |
| **src/evaluator.py** | `dock_molecule()` 完成时调用 `logger.log_event("docking_progress", ...)` |
| **molcraft_agent/tools.py** | `run_synthesis_planning` 中记录 `molecule` 事件<br>`report_iteration` 中记录 `hypothesis_validation` 事件 |

## 4. 关键设计决策

### 4.1 JSON Lines vs 单次写入 JSON
- **选择 JSON Lines**：运行中可流式追加，任何时刻文件有效
- **优势**：运行中断不会产生无效 JSON，支持实时 tail 监控
- **平台兼容**：已验证平台接受 JSON Lines 格式

### 4.2 终端 vs 文件输出分离
- **终端**：保持人类可读的文本格式，不影响实时监控
- **文件**：纯 JSON Lines 格式，专为平台评分优化
- **优势**：不牺牲开发体验，同时满足评分要求

### 4.3 自动打包时机
- **位置**：main.py 的 finally 块
- **时机**：无论成功失败，只要 Agent 执行结束就打包
- **内容**：result.csv + result.log

## 5. 向后兼容

- 终端输出保持原有风格，开发者监控体验不变
- `docs/iteration_log.jsonl` 和 `experiment_round_X.md` 保留不变
- 原有工具链（convert_log_to_json.py）可保留作为备用

## 6. 验收标准

- [ ] `result.log` 每一行都是有效的 JSON 对象
- [ ] 运行过程中可实时查看日志（tail -f result.log）
- [ ] 终端输出保持人类可读
- [ ] Agent 结束后自动生成 result.zip
- [ ] zip 包内包含 result.csv 和 result.log
- [ ] 平台评分系统可正常解析并评分

## 7. 实施计划（下一步）

调用 writing-plans 技能生成详细实施计划，包括：
1. 具体代码修改
2. 测试方案
3. 回滚方案
