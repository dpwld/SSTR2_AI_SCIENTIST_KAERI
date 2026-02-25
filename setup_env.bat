@echo off
chcp 65001 >nul 2>&1
echo.
echo ============================================================
echo   SSTR2 AI Scientist - Windows 환경 설정
echo   (venv + torch + PeptiVerse + ESM-2 + 전체 패키지)
echo ============================================================
echo.

REM --- Python 확인 ---
python --version >nul 2>&1
if errorlevel 1 (
    echo [FAIL] Python이 설치되어 있지 않습니다.
    echo        https://www.python.org/downloads/ 에서 설치하세요.
    pause
    exit /b 1
)

echo [OK] Python 감지:
python --version

REM --- 설정 스크립트 실행 ---
echo.
echo 환경 설정을 시작합니다...
echo   - 가상환경 생성
echo   - PyTorch 설치 (GPU 자동 감지)
echo   - 전체 패키지 설치
echo   - PeptiVerse / ESM-2 모델 다운로드
echo.

python setup_env.py %*

if errorlevel 1 (
    echo.
    echo [WARN] 일부 설치에 실패했습니다. 위 로그를 확인하세요.
) else (
    echo.
    echo [OK] 환경 설정 완료!
)

echo.
echo ============================================================
echo   다음 단계:
echo   1. .venv\Scripts\activate.bat    (가상환경 활성화)
echo   2. ollama serve                  (별도 터미널에서)
echo   3. python run_pipeline.py --test-mode  (파이프라인 실행)
echo   4. streamlit run dashboard.py    (대시보드)
echo ============================================================
echo.
pause
