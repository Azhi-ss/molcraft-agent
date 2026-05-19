# MolCraft Agent

MolCraft Agent 是一个基于 Kimi Agent SDK 的自主科研智能体，用于靶向药物小分子设计与合成路线规划。它遵循文献解析→瓶颈诊断→代码演进→实验验证四阶段迭代循环，自主完成药物发现流程。

## Language

**Session**:
一次 Agent 运行的完整对话上下文，包含 program.md、所有 user/assistant/tool 消息。持久化存储在 `~/.kimi/sessions/` 下，可以跨运行恢复。
_Avoid_: 会话、对话

**Resume**:
从上次中断的 Session 接续执行，不重建对话上下文。通过 `--resume` CLI 参数触发。
_Avoid_: 恢复、重连

**Program**:
定义 Agent 行为规范和工作流的指令文档（`program.md`），在 Session 创建时作为首条 user 消息注入。
_Avoid_: 指令书、任务定义

**Autoresearch Loop**:
Agent 的自主科研循环：文献解析 → 瓶颈诊断 → 代码演进 → 实验验证。每完成一轮称为一次 Iteration。
_Avoid_: 工作流、管线

**Auto-approval (YOLO)**:
Agent 自动批准所有操作请求，无需人工确认。
_Avoid_: 无人值守

**Iteration**:
Autoresearch Loop 的一轮完整执行，包含全部四个阶段。Agent 通过 `report_iteration` 工具记录。
_Avoid_: 轮次

**Stage**:
Iteration 内部的四个阶段之一：诊断、代码演进、实验验证、复盘。通过 `begin_stage`/`end_stage` 工具标记。
_Avoid_: 阶段

**Checkpoint**:
kimi-cli 在 context 中自动插入的状态标记点，用于 context 压缩和回滚。写入 `context.jsonl` 的 `_checkpoint` 消息。
_Avoid_: 保存点

**Context Compaction**:
kimi-cli 内置的自动 context 压缩机制。当 token 数量接近 `max_context_size` 时，用 LLM 总结压缩历史消息，保持对话窗口不溢出。
_Avoid_: 压缩

**Auto-detection**:
`--resume` 通过 kimi-cli 底层 `continue_()` 自动定位最近一次 session，无需手动指定 session_id。

## Relationships

- 一个 **Program** 驱动一个 **Session** 的完整生命周期
- 一个 **Session** 可以被多次 **Resume**，每次恢复追加一条 user 消息
- 一个 **Autoresearch Loop** 包含多个 **Iterations**
- 一个 **Iteration** 包含 4 个 **Stages**
- **Context Compaction** 在 token 数接近上限时自动触发，不影响 **Resume**
- **Checkpoints** 在每条 user/assistant/tool 消息插入后自动创建

## Example dialogue

> **Dev:** "Agent 在迭代 2 的代码演进阶段被 90 分钟超时杀了。我可以用 `--resume` 让它接上吗？"
>
> **Domain expert:** "可以。`--resume` 会找到上次的 **Session**，加载全部对话历史，发一条接续指令。Agent 会读取 context 判断自己执行到哪一步，跳过已完成阶段继续迭代。注意：如果 context 过大，**Context Compaction** 会压缩早期消息再恢复。"

> **Dev:** "如果 Agent 已经跑完了 3 轮迭代，我误操作 `--resume` 会怎样？"
>
> **Domain expert:** "Agent 重新加载 context 发现 '已完成所有迭代'，会进入最终输出阶段。不会造成破坏，只是重复输出而已。目前不检测完成状态。"

## Flagged ambiguities

- "恢复" 曾同时用于 "从代码错误恢复" 和 "从 session 断点恢复" → 已解决：后者用英文 **Resume** 作为规范术语
- "session" 曾混用为 "kimi-cli 的底层 session" 和 "一次 agent 运行" → 已解决：统一指 **Session**