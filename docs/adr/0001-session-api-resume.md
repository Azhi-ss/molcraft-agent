# ADR-0001: 降级至低层 Session API 以支持断点续传

需要用 `--resume` 实现中断后接续执行。高层 `prompt()` API 不暴露 session 生命周期管理，无法实现 Session.resume()。将主循环从 `prompt()` 迁移至 `Session.create()`/`Session.resume()` + `session.prompt()`，统一两套路径。

## Status

accepted

## Context

当前 main.py 使用 `kimi_agent_sdk.prompt()` 高层 API，内部自动创建 Session、运行、关闭。每次执行完后 session 数据持久化在 `~/.kimi/sessions/` 下，但：
1. `prompt()` 不返回 session_id，外部无法跟踪当前 session
2. kimi-cli 底层 `WorkDirMeta.last_session_id` 未被更新，`continue_()` 无法定位最新 session
3. 无 API 入口传入现有 session 进行恢复

## Considered Options

- **方案 A（当前方案 + --resume 走低层 API）**：正常启动继续用 `prompt()`，`--resume` 单独走 `Session.resume()` + `session.prompt()`。改动最小，但两套代码路径维护成本略高。
- **方案 B（全部降级为低层 Session API）**：无论是新建还是恢复，都统一走 `Session.create()` / `Session.resume()` + `session.prompt()`。消除分支，但代码从一行变成 ~30 行。

选择了方案 B，因为 `Session.create()` 的参数签名与 `prompt()` 几乎完全一致，降级成本很低，统一路径避免了 "为什么新建和恢复用不同代码" 的困惑。

## Consequences

- 新建 session 时，`last_session_id` 会被 kimi-cli 底层自动更新（低层 API 保留此行为），`--resume` 可通过 `continue_()` 定位
- `session.prompt()` 返回 WireMessage 而非聚合后的 Message，需要额外处理 ApprovalRequest（已由 yolo=True 自动处理）
- 恢复时发送 ~400 字符的接续指令而非完整 program.md（已在 context 中）
- 可以通过 `session.id` 属性暴露当前 session_id 给日志