#!/bin/bash
# MolCraft Agent + PocketXMol 扩散模型启动脚本
# 自动管理 GPU 服务器 SSH 隧道

set -e

GPU_SERVER="gpu-server"
TUNNEL_PORT="${DIFFUSION_TUNNEL_PORT:-8001}"
GPU_PORT="${GPU_PORT:-8000}"
PXM_SERVER_LOG="/tmp/pxm_server.log"

echo "=== MolCraft Agent + PocketXMol ==="

# 1. 检查 SSH 配置
echo "[1/4] 检查 SSH 连接..."
if ! ssh -o ConnectTimeout=5 "$GPU_SERVER" 'echo ok' > /dev/null 2>&1; then
    echo "警告: 无法 SSH 连接到 $GPU_SERVER, 扩散模型功能将不可用"
else
    echo "   SSH 连接正常"

    # 2. 确保 GPU 服务器上的 PocketXMol 服务正在运行
    echo "[2/4] 检查 PocketXMol 服务..."
    if ssh "$GPU_SERVER" 'pgrep -f "server/main.py" > /dev/null 2>&1'; then
        echo "   PocketXMol 服务已在运行"
    else
        echo "   正在启动 PocketXMol 服务..."
        ssh "$GPU_SERVER" "cd /root/PocketXMol && nohup python server/main.py --pxm-dir /root/PocketXMol --port $GPU_PORT --device cuda:0 > $PXM_SERVER_LOG 2>&1 &"
        echo "   等待服务启动..."
        for i in {1..10}; do
            if ssh "$GPU_SERVER" "curl -s http://localhost:$GPU_PORT/health | grep -q ok"; then
                echo "   PocketXMol 服务已启动"
                break
            fi
            sleep 2
        done
    fi

    # 3. 建立 SSH 隧道
    echo "[3/4] 建立 SSH 隧道 localhost:$TUNNEL_PORT -> $GPU_SERVER:$GPU_PORT..."
    pkill -f "ssh.*$TUNNEL_PORT.*localhost:$GPU_PORT.*$GPU_SERVER" 2>/dev/null || true
    ssh -fNL "$TUNNEL_PORT:localhost:$GPU_PORT" "$GPU_SERVER"
    sleep 1
    if curl -s "http://localhost:$TUNNEL_PORT/health" | grep -q ok; then
        echo "   隧道建立成功, 扩散模型可用"
    else
        echo "   隧道建立失败, 扩散模型功能将不可用"
    fi
fi

# 4. 启动 Agent
echo "[4/4] 启动 MolCraft Agent..."
echo ""

# 捕获退出信号，清理隧道
cleanup() {
    echo ""
    echo "=== 清理 ==="
    pkill -f "ssh.*$TUNNEL_PORT.*localhost:$GPU_PORT.*$GPU_SERVER" 2>/dev/null || true
    echo "SSH 隧道已关闭"
}
trap cleanup EXIT

# 使用 source 加载 .venv 环境
source .venv/bin/activate
python main.py "$@"
