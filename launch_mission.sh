#!/bin/bash
# ==========================================================================
#  SSTR2 AI Scientist — Tmux Mission Control Launcher
# ==========================================================================
#
#  Windows Git Bash + Tmux 환경에서 실행.
#  터미널을 3개 Pane으로 분할:
#    Pane 0 (좌측 상단): 메인 파이프라인 실행 (Module 00~03)
#    Pane 1 (우측 상단): Streamlit 대시보드 서버
#    Pane 2 (하단):       Boardroom 로그 실시간 모니터링
#
#  [사용법]
#    bash launch_mission.sh
#
#  [사전 요구사항]
#    - tmux 설치됨 (Git Bash: pacman -S tmux, 또는 MSYS2)
#    - Python 가상환경 활성화 상태 (conda activate sstr2)
#    - pip install -r requirements.txt 완료
# ==========================================================================

SESSION_NAME="sstr2-mission"
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"

# 기존 세션이 있으면 종료
tmux has-session -t "$SESSION_NAME" 2>/dev/null && tmux kill-session -t "$SESSION_NAME"

echo "=========================================="
echo "  SSTR2 AI Scientist - Mission Control"
echo "  Launching tmux session: $SESSION_NAME"
echo "=========================================="

# ---- 새 세션 생성 (Pane 0: 메인 파이프라인) ----
tmux new-session -d -s "$SESSION_NAME" -c "$PROJECT_DIR" -x 200 -y 50

# Pane 0: 메인 파이프라인 순차 실행
tmux send-keys -t "$SESSION_NAME" "echo '=== SSTR2 AI Scientist Pipeline ===' && \
echo 'Starting Module 00: Peptide Generation...' && \
PYTHONIOENCODING=utf-8 python 00_peptide_gen.py && \
echo '' && echo 'Starting Module 01: PK Filter...' && \
PYTHONIOENCODING=utf-8 python 01_fast_pk_filter.py && \
echo '' && echo 'Starting Module 02: Surrogate Model...' && \
PYTHONIOENCODING=utf-8 python 02_train_surrogate.py && \
echo '' && echo 'Starting Module 03: AI Boardroom...' && \
PYTHONIOENCODING=utf-8 python 03_ai_boardroom.py && \
echo '' && echo '=== ALL MODULES COMPLETE ===' || \
echo '=== PIPELINE ERROR ===' " Enter

# ---- Pane 1 (우측): Streamlit 대시보드 ----
tmux split-window -h -t "$SESSION_NAME" -c "$PROJECT_DIR"
tmux send-keys -t "$SESSION_NAME" "echo 'Starting Streamlit Dashboard...' && \
sleep 2 && \
PYTHONIOENCODING=utf-8 streamlit run dashboard.py --server.port 8501 --server.headless true" Enter

# ---- Pane 2 (하단): 로그 모니터링 ----
tmux split-window -v -t "$SESSION_NAME:0.0" -c "$PROJECT_DIR"
tmux send-keys -t "$SESSION_NAME" "echo 'Waiting for Boardroom log...' && \
while [ ! -f output/Boardroom_Transcript.log ]; do sleep 2; echo 'Waiting...'; done && \
echo '=== Log streaming started ===' && \
tail -f output/Boardroom_Transcript.log" Enter

# ---- 레이아웃 정리: 좌측 2/3, 우측 1/3 ----
tmux select-layout -t "$SESSION_NAME" main-vertical

# Pane 0(파이프라인)에 포커스
tmux select-pane -t "$SESSION_NAME:0.0"

# ---- 세션 접속 ----
echo ""
echo "  Tmux session '$SESSION_NAME' launched!"
echo "  Pane 0 (top-left):  Pipeline execution"
echo "  Pane 1 (right):     Streamlit dashboard (http://localhost:8501)"
echo "  Pane 2 (bottom-left): Log tailing"
echo ""
echo "  Attaching to session..."
echo "  (Ctrl+B, D to detach | Ctrl+B, arrow to switch panes)"
echo ""

tmux attach-session -t "$SESSION_NAME"
