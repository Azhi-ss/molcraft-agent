# result.log JSON 格式修复与事件补齐设计

**日期**: 2026-05-18
**版本**: 1.0
**状态**: 待审核

## 1. 背景

### 1.1 现存问题

- `result.log` 中有两个不同形状的 `metrics` 事件：`pipeline.py` 产出 `{molecule_count, trivial_count, avg_binding_energy, min_binding_energy}`，`main.py` 产出 `{sample_count, trivial_ratio}`，字段名不一致
- 原始设计文档定义的 `stage`、`molecule` 事件从未实现
- 日志逻辑散布在 `main.py` 和 `pipeline.py` 的业务代码中，难以单独测试
- 没有单元测试覆盖日志输出格式

### 1.2 修复范围

1. 统一 `metrics` 事件字段形状
2. 补齐 `stage` 事件（Agent 层面 4 个阶段：诊断、代码演进、实验验证、复盘）
3. 补齐 `molecule` 事件（每分子详细数据，含 docking_std）
4. 提取独立 `EventLogger` 模块，分离日志逻辑与业务逻辑
5. 用 TDD 方式驱动实现

## 2. 架构设计

### 2.1 新增文件

| 文件 | 职责 |
|------|------|
| `src/event_schema.py` | 事件类型 dataclass 定义 + 字段校验 |
| `src/event_logger.py` | EventLogger 类，封装 JSONL 写入 + 终端输出 |
| `tests/test_event_schema.py` | Schema 验证测试 |
| `tests/test_event_logger.py` | EventLogger 单元测试 |

### 2.2 修改文件

| 文件 | 改动 |
|------|------|
| `main.py` | `StructuredLogger` → `EventLogger`；删除 `write_result_log` 合并逻辑 |
| `tools/pipeline.py` | 删除 `write_result_log()`；`run_evolutionary_pipeline()` 接受 `EventLogger` 参数 |

### 2.3 模块依赖关系

```
main.py
  ↓ 实例化
EventLogger ──→ event_schema.py (dataclass 定义)
  ↓ 传入
pipeline.py ──→ 调用 logger.log_molecule(), logger.log_docking_progress()
```

`event_logger.py` 和 `event_schema.py` 是独立模块，不依赖任何业务代码。

## 3. 事件 Schema 定义

### StageEvent

```python
@dataclass
class StageEvent:
    name: str              # "诊断" | "代码演进" | "实验验证" | "复盘"
    status: str            # "begun" | "completed" | "failed"
    duration_seconds: float | None = None
```

### MoleculeEvent

```python
@dataclass
class MoleculeEvent:
    mol_smiles: str                       # 必填
    binding_energy: float                 # 必填
    composite_score: float                # 必填
    syn_steps: int                        # 必填
    sa_score: float | None                # 必填（可为 None 如果计算失败）
    qed: float                            # 必填
    rings: int                            # 必填
    trivial: bool = False                 # 可选，可从 syn_steps == 0 推导
    route_quality: float | None = None    # 可选，规划失败时为 None
    docking_std: float | None = None      # 可选，共识对接标准差
```

### MetricsEvent（统一形状）

```python
@dataclass
class MetricsEvent:
    molecule_count: int
    non_trivial_count: int
    trivial_count: int                    # 可从 molecule_count - non_trivial_count 推导但保留
    avg_binding_energy: float | None
    min_binding_energy: float | None
    avg_syn_steps: float | None
    docking_success_rate: float
```

### StartEvent / EndEvent

保持现有字段不变，只做 dataclass 封装。

### DockingProgressEvent

保持现有字段不变。

### HypothesisValidationEvent

保持现有字段不变。

## 4. EventLogger API

```python
class EventLogger:
    def __init__(self, log_path: Path)
    def commit(self)
    def close(self)
    def write(self, message: str)       # 终端输出；文件写 stdout 事件

    # 结构化事件
    def log_start(self, round, max_minutes, max_iterations, program_file)
    def log_end(self, status)
    def log_stage(self, name: str, status: str, duration_seconds: float | None = None)
    def log_docking_progress(self, current, total, success_rate)
    def log_molecule(self, event: MoleculeEvent)
    def log_metrics(self, event: MetricsEvent)
    def log_hypothesis_validation(self, hypothesis_id, success, conclusion, changes=None)
```

### Stage 事件的触发方式

Agent 工作流是 LLM 驱动的，阶段切换发生在 LLM 推理和工具调用中，无法在 `main.py` 用 context manager 硬编码包裹。

因此 stage 事件通过 **Agent 工具调用** 触发：
- `tools.py` 中新增 `begin_stage(stage_name)` / `end_stage()` 工具函数
- `program.md` 的各阶段开头/末尾指示 Agent 调用对应工具
- `end_stage()` 自动计算耗时并写入 `stage` 事件

```python
# tools.py
def begin_stage(name: str):
    """Agent 调用：标记一个阶段的开始。"""
    _stage_name = name
    _stage_start = time.monotonic()
    sys.stdout.log_event("stage", name=name, status="begun")
    return f"阶段「{name}」开始"

def end_stage():
    """Agent 调用：标记当前阶段的结束，自动计算耗时。"""
    elapsed = time.monotonic() - _stage_start
    sys.stdout.log_event("stage", name=_stage_name, status="completed",
                         duration_seconds=round(elapsed, 1))
    return f"阶段「{_stage_name}」完成，耗时 {elapsed:.1f}s"
```

## 5. JSONL 输出示例

```
{"type": "start", "timestamp": "2026-05-18T21:19:51", "round": 10, "max_minutes": 120, "max_iterations": 3, "program_file": "program.md"}
{"type": "stage", "timestamp": "2026-05-18T21:19:52", "name": "诊断", "status": "begun"}
{"type": "stdout", "timestamp": "2026-05-18T21:19:52", "content": "## 诊断瓶颈..."}
...
{"type": "stage", "timestamp": "2026-05-18T21:25:10", "name": "诊断", "status": "completed", "duration_seconds": 318.0}
{"type": "stage", "timestamp": "2026-05-18T21:25:11", "name": "代码演进", "status": "begun"}
...
{"type": "stage", "timestamp": "2026-05-18T21:35:00", "name": "代码演进", "status": "completed", "duration_seconds": 589.0}
{"type": "stage", "timestamp": "2026-05-18T21:35:01", "name": "实验验证", "status": "begun"}
{"type": "docking_progress", "timestamp": "...", "current": 5, "total": 50, "success_rate": 1.0}
...
{"type": "molecule", "timestamp": "...", "mol_smiles": "Cc1cc(C)c2ccccc2c1O", "binding_energy": -8.098, "composite_score": 0.97, "syn_steps": 2, "sa_score": 3.2, "qed": 0.68, "rings": 3, "trivial": false, "route_quality": 0.85, "docking_std": 0.05}
...
{"type": "stage", "timestamp": "2026-05-18T22:04:25", "name": "实验验证", "status": "completed", "duration_seconds": 1764.0}
{"type": "metrics", "timestamp": "...", "molecule_count": 10, "non_trivial_count": 10, "trivial_count": 0, "avg_binding_energy": -8.314, "min_binding_energy": -9.074, "avg_syn_steps": 1.5, "docking_success_rate": 0.95}
{"type": "stage", "timestamp": "2026-05-18T22:04:26", "name": "复盘", "status": "begun"}
...
{"type": "stage", "timestamp": "2026-05-18T22:09:19", "name": "复盘", "status": "completed", "duration_seconds": 293.0}
{"type": "end", "timestamp": "2026-05-18T22:09:19", "status": "success"}
```

## 6. TDD 测试计划

### 第 1 层：Schema 测试 (`tests/test_event_schema.py`)

- 缺必填字段抛 ValueError
- 可选字段默认值正确
- `to_dict()` 序列化不丢字段
- 一致性校验：`trivial_count + non_trivial_count == molecule_count`

### 第 2 层：EventLogger 单元测试 (`tests/test_event_logger.py`)

- 各类事件写入有效 JSONL
- `log_stage("诊断", "completed", duration_seconds=318)` 输出正确
- commit 原子替换
- close 不 commit 清理 tmp
- NaN 拒绝
- 终端输出纯文本

### 第 2b 层：tools.py begin_stage/end_stage 测试

- `begin_stage("诊断")` 写入 `{"status": "begun"}` 事件
- `end_stage()` 写入 `{"status": "completed", "duration_seconds": ...}` 事件
- `begin_stage` 后未 `end_stage` 时第二次调用抛错（防嵌套）

### 第 3 层：集成测试

- 运行后 result.log 每行 valid JSON
- 所有事件类型齐全
- metrics 字段与 CSV 一致

## 7. 向后兼容

- `docs/iteration_log.jsonl` 和 `experiment_round_X.md` 不受影响
- 终端输出保持人类可读
- `MOLCRAFT_LOG_PATH` 环境变量机制移除，改为直接传递 EventLogger 实例

## 8. 验收标准

- [ ] `result.log` 中只有一个 `metrics` 事件形状
- [ ] `result.log` 包含 `stage` 事件（Agent 4 阶段）
- [ ] `result.log` 包含 `molecule` 事件（每分子 7 必选字段 + 可选字段）
- [ ] `molecule` 事件包含 `docking_std`
- [ ] 所有测试通过
- [ ] 终端输出保持人类可读
- [ ] `result.zip` 正常生成
