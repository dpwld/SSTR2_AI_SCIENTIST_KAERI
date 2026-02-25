@echo off
chcp 65001 >nul 2>&1
echo.
echo ============================================================
echo   SSTR2 AI Scientist - Pipeline Runner (Windows)
echo ============================================================
echo.

REM --- 가상환경 자동 활성화 ---
if exist ".venv\Scripts\activate.bat" (
    echo [INFO] 가상환경 활성화 중...
    call .venv\Scripts\activate.bat
) else (
    echo [WARN] 가상환경 없음. 시스템 Python을 사용합니다.
    echo        먼저 setup_env.bat 를 실행하세요.
    echo.
)

REM --- 파이프라인 실행 ---
set PYTHONIOENCODING=utf-8
python run_pipeline.py --test-mode %*

echo.
pause
