#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  SSTR2 AI Scientist — 원클릭 환경 설정 스크립트  (setup_env.py)
==========================================================================

  이 스크립트는 프로젝트 실행에 필요한 모든 의존성을 자동으로 설치하고
  환경을 검증합니다. Python 3.10+ 환경에서 실행하세요.

  [사용법]
    python setup_env.py               # 기본 설치 (venv + pip + torch CPU + 모델 다운로드)
    python setup_env.py --no-venv     # 가상환경 없이 현재 환경에 직접 설치
    python setup_env.py --cpu-only    # torch GPU 감지 생략, CPU 버전만 설치
    python setup_env.py --full        # 전체 설치 + Ollama 모델 pull

  [설치 순서]
    1. Python 버전 확인 (>=3.10)
    2. 가상환경(venv) 생성 및 활성화
    3. pip 업그레이드
    4. torch 설치 (NVIDIA GPU 자동 감지 → CUDA / CPU 분기)
    5. requirements.txt 나머지 패키지 설치
    6. PeptiVerse + ESM-2 모델 가중치 사전 다운로드
    7. Ollama 설치 및 모델 확인
    8. 전체 환경 검증 리포트 출력
==========================================================================
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

# ─── 색상 출력 헬퍼 ─────────────────────────────────────────────────────
IS_WIN = platform.system() == "Windows"

def _c(text, color):
    """ANSI 색상 코드로 텍스트 감싸기 (Windows 터미널 호환)."""
    colors = {"green": "32", "red": "31", "yellow": "33", "cyan": "36", "bold": "1"}
    if IS_WIN:
        os.system("")  # Windows ANSI escape 활성화
    code = colors.get(color, "0")
    return f"\033[{code}m{text}\033[0m"

def ok(msg):    print(f"  {_c('[OK]', 'green')}   {msg}")
def warn(msg):  print(f"  {_c('[WARN]', 'yellow')} {msg}")
def fail(msg):  print(f"  {_c('[FAIL]', 'red')} {msg}")
def info(msg):  print(f"  {_c('[INFO]', 'cyan')} {msg}")

# ─── 프로젝트 경로 ──────────────────────────────────────────────────────
PROJECT_DIR = Path(__file__).resolve().parent
VENV_DIR = PROJECT_DIR / ".venv"
REQ_FILE = PROJECT_DIR / "requirements.txt"

# ─── Windows/Linux/Mac 분기 ─────────────────────────────────────────────
if IS_WIN:
    VENV_PYTHON = VENV_DIR / "Scripts" / "python.exe"
    VENV_PIP = VENV_DIR / "Scripts" / "pip.exe"
    VENV_ACTIVATE = VENV_DIR / "Scripts" / "activate.bat"
    ACTIVATE_CMD = str(VENV_ACTIVATE)
else:
    VENV_PYTHON = VENV_DIR / "bin" / "python"
    VENV_PIP = VENV_DIR / "bin" / "pip"
    VENV_ACTIVATE = VENV_DIR / "bin" / "activate"
    ACTIVATE_CMD = f"source {VENV_ACTIVATE}"


# ==========================================================================
#  Step 1: Python 버전 확인
# ==========================================================================

def check_python_version():
    """Python 3.10 이상인지 확인."""
    print("\n" + "=" * 60)
    print("  Step 1: Python 버전 확인")
    print("=" * 60)

    v = sys.version_info
    ver_str = f"{v.major}.{v.minor}.{v.micro}"

    if v.major < 3 or (v.major == 3 and v.minor < 10):
        fail(f"Python {ver_str} 감지. Python 3.10 이상이 필요합니다.")
        fail("https://www.python.org/downloads/ 에서 최신 버전을 설치하세요.")
        sys.exit(1)

    ok(f"Python {ver_str} ({platform.python_implementation()})")
    ok(f"실행 경로: {sys.executable}")
    return True


# ==========================================================================
#  Step 2: 가상환경 생성
# ==========================================================================

def setup_venv(skip_venv=False):
    """가상환경(.venv)을 생성한다."""
    print("\n" + "=" * 60)
    print("  Step 2: 가상환경(venv) 설정")
    print("=" * 60)

    if skip_venv:
        info("--no-venv 옵션: 가상환경 생략, 현재 환경에 직접 설치합니다.")
        return sys.executable

    if VENV_DIR.exists():
        info(f"기존 가상환경 감지: {VENV_DIR}")
        if VENV_PYTHON.exists():
            ok(f"가상환경 Python: {VENV_PYTHON}")
            return str(VENV_PYTHON)
        else:
            warn("가상환경이 손상됨. 재생성합니다...")
            shutil.rmtree(VENV_DIR)

    info(f"가상환경 생성 중: {VENV_DIR}")
    subprocess.run([sys.executable, "-m", "venv", str(VENV_DIR)], check=True)

    if not VENV_PYTHON.exists():
        fail("가상환경 생성 실패!")
        sys.exit(1)

    ok(f"가상환경 생성 완료: {VENV_DIR}")
    ok(f"활성화 명령: {ACTIVATE_CMD}")
    return str(VENV_PYTHON)


# ==========================================================================
#  Step 3: NVIDIA GPU 감지
# ==========================================================================

def detect_nvidia_gpu():
    """NVIDIA GPU + CUDA 사용 가능 여부를 확인한다."""
    try:
        result = subprocess.run(
            ["nvidia-smi"], capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            # CUDA 버전 파싱 시도
            for line in result.stdout.split("\n"):
                if "CUDA Version" in line:
                    info(f"  {line.strip()}")
                    break
            return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return False


# ==========================================================================
#  Step 4: torch 설치 (GPU 자동 감지)
# ==========================================================================

def install_torch(python_exe, cpu_only=False):
    """PyTorch를 NVIDIA GPU 유무에 따라 CUDA/CPU 버전으로 설치한다."""
    print("\n" + "=" * 60)
    print("  Step 3: PyTorch 설치 (PeptiVerse/ESM-2 임베딩 엔진)")
    print("=" * 60)

    # 이미 설치 확인
    check = subprocess.run(
        [python_exe, "-c", "import torch; print(torch.__version__)"],
        capture_output=True, text=True
    )
    if check.returncode == 0:
        ver = check.stdout.strip()
        ok(f"torch 이미 설치됨: v{ver}")
        # CUDA 확인
        cuda_check = subprocess.run(
            [python_exe, "-c",
             "import torch; print('CUDA' if torch.cuda.is_available() else 'CPU')"],
            capture_output=True, text=True
        )
        device = cuda_check.stdout.strip() if cuda_check.returncode == 0 else "?"
        ok(f"torch 디바이스: {device}")
        return True

    # GPU 감지
    has_gpu = False
    if not cpu_only:
        info("NVIDIA GPU 감지 중...")
        has_gpu = detect_nvidia_gpu()

    if has_gpu and not cpu_only:
        ok("NVIDIA GPU 감지 → CUDA 버전 torch 설치")
        info("torch + CUDA 12.1 설치 중 (약 2~5GB 다운로드)...")
        result = subprocess.run(
            [python_exe, "-m", "pip", "install",
             "torch", "--index-url", "https://download.pytorch.org/whl/cu121"],
            timeout=1800
        )
    else:
        if cpu_only:
            info("--cpu-only 옵션 → CPU 버전 torch 설치")
        else:
            info("NVIDIA GPU 미감지 → CPU 버전 torch 설치")
        info("torch (CPU) 설치 중 (약 200MB 다운로드)...")
        result = subprocess.run(
            [python_exe, "-m", "pip", "install", "torch", "--index-url",
             "https://download.pytorch.org/whl/cpu"],
            timeout=1800
        )

    if result.returncode != 0:
        warn("torch 설치 실패. pip install torch 수동 재시도 중...")
        subprocess.run(
            [python_exe, "-m", "pip", "install", "torch"],
            timeout=1800
        )

    # transformers 설치
    info("transformers 설치 중...")
    subprocess.run(
        [python_exe, "-m", "pip", "install", "transformers>=4.30.0"],
        capture_output=True, timeout=600
    )

    # 설치 확인
    check2 = subprocess.run(
        [python_exe, "-c", "import torch; print(torch.__version__)"],
        capture_output=True, text=True
    )
    if check2.returncode == 0:
        ok(f"torch 설치 완료: v{check2.stdout.strip()}")
        return True
    else:
        fail("torch 설치 실패! 수동 설치가 필요합니다.")
        info("  pip install torch")
        info("  pip install transformers")
        return False


# ==========================================================================
#  Step 5: 나머지 패키지 설치
# ==========================================================================

def install_requirements(python_exe):
    """requirements.txt 패키지를 설치한다 (torch는 이미 설치됨)."""
    print("\n" + "=" * 60)
    print("  Step 4: 나머지 패키지 설치")
    print("=" * 60)

    # pip 업그레이드
    info("pip 업그레이드 중...")
    subprocess.run(
        [python_exe, "-m", "pip", "install", "--upgrade", "pip"],
        capture_output=True
    )

    # requirements.txt 설치 (torch는 이미 있으므로 skip됨)
    info(f"requirements.txt 패키지 설치 중...")
    result = subprocess.run(
        [python_exe, "-m", "pip", "install", "-r", str(REQ_FILE)],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        warn("일부 패키지 설치 실패. 개별 재시도 중...")
        print(result.stderr[-500:] if result.stderr else "")

        # 실패하기 쉬운 패키지 개별 설치
        for pkg in ["rdkit-pypi", "langgraph", "langchain-ollama", "langchain-core"]:
            info(f"  {pkg} 개별 설치 시도...")
            subprocess.run(
                [python_exe, "-m", "pip", "install", pkg],
                capture_output=True
            )
    else:
        ok("requirements.txt 패키지 설치 완료")


# ==========================================================================
#  Step 6: PLM 모델 가중치 사전 다운로드
# ==========================================================================

def predownload_models(python_exe):
    """PeptiVerse와 ESM-2 모델 가중치를 HuggingFace에서 사전 다운로드한다.

    첫 실행 시 다운로드 지연을 방지하기 위해 미리 캐시에 저장.
    """
    print("\n" + "=" * 60)
    print("  Step 5: PLM 모델 가중치 사전 다운로드")
    print("=" * 60)

    # torch/transformers 설치 확인
    check = subprocess.run(
        [python_exe, "-c", "import torch; import transformers"],
        capture_output=True, text=True
    )
    if check.returncode != 0:
        warn("torch/transformers 미설치 → 모델 다운로드 건너뜀")
        info("수동 특성 공학 Fallback으로 파이프라인 실행 가능")
        return

    # 1순위: PeptiVerse
    info("1순위: PeptiVerse (ChatterjeeLab/PeptiVerse) 다운로드 중...")
    peptiverse_script = (
        "from transformers import AutoTokenizer, AutoModel; "
        "t = AutoTokenizer.from_pretrained('ChatterjeeLab/PeptiVerse', trust_remote_code=True); "
        "m = AutoModel.from_pretrained('ChatterjeeLab/PeptiVerse', trust_remote_code=True); "
        "print(f'hidden_size={m.config.hidden_size}')"
    )
    result = subprocess.run(
        [python_exe, "-c", peptiverse_script],
        capture_output=True, text=True, timeout=600
    )
    if result.returncode == 0:
        ok(f"PeptiVerse 다운로드 완료 ({result.stdout.strip()})")
    else:
        warn("PeptiVerse 다운로드 실패 (네트워크 또는 호환성 문제)")
        if result.stderr:
            # 핵심 에러만 출력
            for line in result.stderr.strip().split("\n")[-3:]:
                info(f"  {line.strip()}")

    # 2순위: ESM-2 (Fallback 모델)
    info("2순위: ESM-2 (facebook/esm2_t6_8M_UR50D) 다운로드 중...")
    esm2_script = (
        "from transformers import AutoTokenizer, AutoModel; "
        "t = AutoTokenizer.from_pretrained('facebook/esm2_t6_8M_UR50D'); "
        "m = AutoModel.from_pretrained('facebook/esm2_t6_8M_UR50D'); "
        "print(f'hidden_size={m.config.hidden_size}')"
    )
    result = subprocess.run(
        [python_exe, "-c", esm2_script],
        capture_output=True, text=True, timeout=600
    )
    if result.returncode == 0:
        ok(f"ESM-2 다운로드 완료 ({result.stdout.strip()})")
    else:
        warn("ESM-2 다운로드 실패")

    info("모델 캐시 위치: ~/.cache/huggingface/hub/")


# ==========================================================================
#  Step 7: Ollama 확인
# ==========================================================================

def check_ollama(pull_model=False):
    """Ollama 설치 및 모델 확인."""
    print("\n" + "=" * 60)
    print("  Step 6: Ollama (로컬 LLM) 확인")
    print("=" * 60)

    # Ollama 설치 확인
    ollama_path = shutil.which("ollama")
    if ollama_path is None:
        warn("Ollama가 설치되어 있지 않습니다.")
        warn("Module 03 (AI 드림팀 원탁회의)에 필요합니다.")
        info("설치 방법:")
        if IS_WIN:
            info("  https://ollama.com/download/windows")
            info("  또는: winget install Ollama.Ollama")
        elif platform.system() == "Darwin":
            info("  https://ollama.com/download/mac")
            info("  또는: brew install ollama")
        else:
            info("  curl -fsSL https://ollama.com/install.sh | sh")
        return False

    ok(f"Ollama 감지: {ollama_path}")

    # Ollama 서버 상태 확인
    try:
        result = subprocess.run(
            ["ollama", "list"], capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            ok("Ollama 서버 응답 확인")
            models = result.stdout.strip()
            if models:
                info("설치된 모델:")
                for line in models.split("\n")[:10]:
                    info(f"  {line}")
            else:
                warn("설치된 모델이 없습니다.")
        else:
            warn("Ollama 서버가 실행 중이 아닙니다.")
            info("실행 방법: ollama serve (별도 터미널에서)")
    except (FileNotFoundError, subprocess.TimeoutExpired):
        warn("Ollama 서버 연결 실패")
        info("실행 방법: ollama serve (별도 터미널에서)")

    # 모델 pull (--full 옵션)
    if pull_model:
        model_name = "llama3.1:8b"
        info(f"모델 다운로드 중: {model_name} (약 4.7GB, 시간이 걸릴 수 있습니다)...")
        result = subprocess.run(
            ["ollama", "pull", model_name],
            timeout=3600  # 최대 1시간
        )
        if result.returncode == 0:
            ok(f"모델 다운로드 완료: {model_name}")
        else:
            warn(f"모델 다운로드 실패. 수동 실행: ollama pull {model_name}")

    return True


# ==========================================================================
#  Step 8: AutoDock Vina 확인 (선택)
# ==========================================================================

def check_vina():
    """AutoDock Vina 설치 확인 (Production 도킹용, 선택사항)."""
    print("\n" + "=" * 60)
    print("  Step 7: AutoDock Vina (선택 — Production 도킹용)")
    print("=" * 60)

    vina_path = shutil.which("vina")
    if vina_path:
        ok(f"AutoDock Vina 감지: {vina_path}")
        return True

    info("AutoDock Vina가 설치되어 있지 않습니다.")
    info("Module 02는 --test-mode로 실행 가능합니다 (더미 도킹 점수).")
    info("실제 도킹을 위한 설치 방법:")
    info("  conda install -c conda-forge autodock-vina")
    info("  또는: https://vina.scripps.edu/ 에서 다운로드")
    return False


# ==========================================================================
#  Step 9: 전체 환경 검증
# ==========================================================================

def verify_environment(python_exe):
    """설치된 패키지와 환경을 종합 검증한다."""
    print("\n" + "=" * 60)
    print("  Step 8: 환경 검증 리포트")
    print("=" * 60)

    # 각 패키지 임포트 테스트
    packages = [
        # (패키지명, 임포트명, 사용처, 필수여부)
        ("numpy", "numpy", "전체 모듈", True),
        ("pandas", "pandas", "전체 모듈", True),
        ("peptides", "peptides", "Module 00, 01", True),
        ("rdkit-pypi", "rdkit", "Module 01 (PPB/Rigidity), 02 (3D 변환)", True),
        ("xgboost", "xgboost", "Module 02 (대리 모델)", True),
        ("scikit-learn", "sklearn", "Module 02 (학습/검증)", True),
        ("torch", "torch", "Module 02 (PeptiVerse/ESM-2 임베딩)", True),
        ("transformers", "transformers", "Module 02 (PeptiVerse/ESM-2 모델)", True),
        ("langgraph", "langgraph", "Module 03 (AI 원탁회의)", True),
        ("langchain-ollama", "langchain_ollama", "Module 03 (Ollama LLM)", True),
        ("langchain-core", "langchain_core", "Module 03 (에이전트)", True),
        ("matplotlib", "matplotlib", "Module 04 (시각화)", True),
        ("seaborn", "seaborn", "Module 04 (통계 시각화)", True),
        ("scipy", "scipy", "Module 04 (Spearman 상관)", True),
        ("streamlit", "streamlit", "Dashboard", True),
        ("plotly", "plotly", "Dashboard (인터랙티브 차트)", True),
    ]

    results = []
    for pkg_name, import_name, usage, required in packages:
        check_cmd = f"import {import_name}; print({import_name}.__version__)"
        result = subprocess.run(
            [python_exe, "-c", check_cmd],
            capture_output=True, text=True
        )

        if result.returncode == 0:
            version = result.stdout.strip()
            ok(f"{pkg_name:<20s} v{version:<12s} [{usage}]")
            results.append((pkg_name, True, version))
        else:
            if required:
                fail(f"{pkg_name:<20s} {'미설치':<12s} [{usage}]")
            else:
                warn(f"{pkg_name:<20s} {'미설치':<12s} [{usage}] (선택사항)")
            results.append((pkg_name, False, None))

    # PeptiVerse 모델 캐시 확인
    print()
    info("PLM 모델 캐시 확인:")
    for model_id, label in [
        ("ChatterjeeLab/PeptiVerse", "PeptiVerse (1순위)"),
        ("facebook/esm2_t6_8M_UR50D", "ESM-2 (2순위)"),
    ]:
        check_script = (
            f"from transformers import AutoConfig; "
            f"c = AutoConfig.from_pretrained('{model_id}', trust_remote_code=True); "
            f"print(c.hidden_size)"
        )
        r = subprocess.run(
            [python_exe, "-c", check_script],
            capture_output=True, text=True, timeout=30
        )
        if r.returncode == 0:
            ok(f"  {label}: 캐시됨 (hidden_size={r.stdout.strip()})")
        else:
            warn(f"  {label}: 미캐시 (첫 실행 시 자동 다운로드)")

    # 요약
    n_ok = sum(1 for _, installed, _ in results if installed)
    n_required = sum(1 for _, _, _, req in packages if req)
    n_required_ok = sum(1 for (_, installed, _), (_, _, _, req)
                        in zip(results, packages) if req and installed)

    print("\n" + "-" * 60)
    print(f"  검증 결과: {n_ok}/{len(packages)} 패키지 설치됨")
    print(f"  필수 패키지: {n_required_ok}/{n_required} 통과")

    if n_required_ok == n_required:
        ok("모든 필수 패키지가 정상 설치되었습니다!")
    else:
        fail("일부 필수 패키지가 누락되었습니다. 위의 [FAIL] 항목을 확인하세요.")

    return n_required_ok == n_required


# ==========================================================================
#  최종 사용법 출력
# ==========================================================================

def print_usage(python_exe, skip_venv=False):
    """설치 완료 후 사용법을 출력한다."""
    print("\n" + "=" * 60)
    print("  설치 완료! 사용법 안내")
    print("=" * 60)

    if not skip_venv and VENV_DIR.exists():
        print(f"""
  [1단계] 가상환경 활성화:
    {ACTIVATE_CMD}

  [2단계] Ollama 서버 시작 (별도 터미널):
    ollama serve

  [3단계] 파이프라인 실행:
    a) 테스트 모드 (도킹 엔진 불필요):
       python run_pipeline.py --test-mode

    b) 개별 모듈 실행:
       python 00_peptide_gen.py
       python 01_fast_pk_filter.py
       python 02_train_surrogate.py --test-mode
       python 03_ai_boardroom.py
       python 04_validate_and_analyze.py

    c) 대시보드:
       streamlit run dashboard.py

  [임베딩 모델 우선순위]:
    1순위: PeptiVerse (ChatterjeeLab/PeptiVerse)
    2순위: ESM-2 (facebook/esm2_t6_8M_UR50D)
    3순위: 수동 특성 공학 (torch 없이도 작동)
""")
    else:
        print(f"""
  파이프라인 실행:
    python run_pipeline.py --test-mode
""")


# ==========================================================================
#  Main
# ==========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="SSTR2 AI Scientist 환경 설정 스크립트"
    )
    parser.add_argument("--no-venv", action="store_true",
                        help="가상환경 없이 현재 환경에 직접 설치")
    parser.add_argument("--cpu-only", action="store_true",
                        help="GPU 감지 생략, CPU 버전 torch만 설치")
    parser.add_argument("--full", action="store_true",
                        help="전체 설치 + Ollama 모델 pull")
    parser.add_argument("--skip-models", action="store_true",
                        help="PLM 모델 사전 다운로드 건너뛰기")
    args = parser.parse_args()

    pull_model = args.full

    print("=" * 60)
    print("  SSTR2 AI Scientist — 환경 설정 마법사")
    print("=" * 60)
    print(f"  프로젝트 경로 : {PROJECT_DIR}")
    print(f"  운영체제       : {platform.system()} {platform.release()}")
    print(f"  Python         : {sys.version}")
    print(f"  torch 모드     : {'CPU only' if args.cpu_only else 'GPU 자동감지'}")
    print(f"  Ollama Pull    : {'예' if pull_model else '아니오'}")

    # Step 1: Python 버전
    check_python_version()

    # Step 2: 가상환경
    python_exe = setup_venv(skip_venv=args.no_venv)

    # Step 3: pip 업그레이드 (venv 내에서)
    subprocess.run(
        [python_exe, "-m", "pip", "install", "--upgrade", "pip"],
        capture_output=True
    )

    # Step 4: torch 선행 설치 (GPU 자동 감지)
    install_torch(python_exe, cpu_only=args.cpu_only)

    # Step 5: 나머지 패키지
    install_requirements(python_exe)

    # Step 6: PLM 모델 사전 다운로드
    if not args.skip_models:
        predownload_models(python_exe)
    else:
        info("--skip-models: 모델 다운로드 건너뜀")

    # Step 7: Ollama
    check_ollama(pull_model=pull_model)

    # Step 8: Vina (정보 제공만)
    check_vina()

    # Step 9: 검증
    all_ok = verify_environment(python_exe)

    # 사용법
    print_usage(python_exe, skip_venv=args.no_venv)

    if all_ok:
        print(_c("\n  환경 설정 완료! 파이프라인을 실행할 준비가 되었습니다.\n", "green"))
    else:
        print(_c("\n  일부 패키지 설치에 실패했습니다. 위 로그를 확인하세요.\n", "yellow"))

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
