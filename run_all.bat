@echo off
chcp 65001 >nul 2>&1
set PYTHONIOENCODING=utf-8
cd /d "%~dp0"

echo ============================================================
echo   SSTR2 AI Scientist - Full Production Pipeline
echo ============================================================
echo.

REM ── Step 1: venv 활성화 ──
echo [Step 1/8] 가상환경 활성화...
if exist ".venv\Scripts\activate.bat" (
    call ".venv\Scripts\activate.bat"
    python --version
) else (
    echo   .venv 없음 - 시스템 Python 사용
)
echo.

REM ── Step 2: 핵심 패키지 개별 설치 확인 ──
echo [Step 2/8] 핵심 패키지 설치 확인...
python install_deps.py
if %ERRORLEVEL% neq 0 (
    echo   [WARNING] 일부 패키지 미설치 - 파이프라인이 fallback으로 진행됩니다.
)
echo.

REM ── Step 3: 수용체 PDBQT 재생성 ──
echo [Step 3/8] 수용체 PDBQT 재생성 (Vina 1.2.x 호환)...
if exist "SSTR2_receptor.pdbqt" del "SSTR2_receptor.pdbqt"
python prepare_receptor.py
if %ERRORLEVEL% neq 0 (
    echo   [ERROR] 수용체 파일 생성 실패!
    pause
    exit /b 1
)
echo.

REM ── Step 4: output 초기화 ──
echo [Step 4/8] output 폴더 초기화...
if exist "output" rmdir /s /q "output"
mkdir "output"
mkdir "output\figures"
echo   완료
echo.

REM ── Step 5: Vina 동작 테스트 ──
echo [Step 5/8] Vina 동작 테스트...
python -c "from config import VINA_EXECUTABLE; import subprocess; r=subprocess.run([VINA_EXECUTABLE,'--version'],capture_output=True,text=True); print(r.stdout.strip() or r.stderr.strip())" 2>&1
echo.

REM ── Step 6: PDBQT 포맷 검증 ──
echo [Step 6/8] 수용체 PDBQT 포맷 검증...
python validate_pdbqt.py
echo.

REM ── Step 7: Ollama 서버 확인 ──
echo [Step 7/8] Ollama 서버 확인...
python -c "import urllib.request; urllib.request.urlopen('http://localhost:11434/api/tags',timeout=3); print('  Ollama OK')" >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo   [경고] Ollama 서버 미감지 - Module 03 건너뛸 수 있음
    echo   시작: ollama serve
) else (
    echo   Ollama 서버 정상
)
echo.

REM ── Step 8: 파이프라인 실행 ──
echo [Step 8/8] 파이프라인 실행 (Production, Module 00~04)...
echo.

python run_pipeline.py
if %ERRORLEVEL% neq 0 (
    echo.
    echo [ERROR] 파이프라인 중 오류 발생!
) else (
    echo.
    echo ============================================================
    echo   파이프라인 완료!
    echo ============================================================
    dir /b "output\*.csv" "output\*.md" "output\*.log" 2>nul
    echo   대시보드: streamlit run dashboard.py
)

echo.
pause
