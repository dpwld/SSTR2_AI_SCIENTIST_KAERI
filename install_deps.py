#!/usr/bin/env python3
"""핵심 패키지 개별 설치 확인 스크립트."""
import subprocess
import sys

if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def check_import(module, from_name=None):
    """모듈 임포트 가능 여부 확인."""
    try:
        if from_name:
            exec(f"from {module} import {from_name}")
        else:
            __import__(module)
        return True
    except ImportError:
        return False


def pip_install(package):
    """pip install 실행."""
    print(f"  설치 중: {package}")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", package],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  [WARNING] {package} 설치 실패: {result.stderr[:200]}")
        return False
    return True


def check_numpy_version():
    """numpy < 2.0 인지 확인."""
    try:
        import numpy
        major = int(numpy.__version__.split(".")[0])
        return major < 2
    except ImportError:
        return False


packages = [
    # (import_name, from_name, pip_package, description)
    ("numpy", None, "numpy>=1.24.0,<2.0.0", "numpy"),
    ("pandas", None, "pandas>=2.0.0", "pandas"),
    ("peptides", None, "peptides>=0.3.0", "peptides"),
    ("rdkit", "Chem", "rdkit", "rdkit"),
    ("xgboost", None, "xgboost>=1.7.0", "xgboost"),
    ("sklearn", None, "scikit-learn>=1.3.0", "scikit-learn"),
    ("torch", None, "torch>=2.0.0", "torch"),
    ("transformers", None, "transformers>=4.30.0", "transformers"),
    ("sentencepiece", None, "sentencepiece>=0.1.99", "sentencepiece"),
    ("meeko", "MoleculePreparation", "meeko", "meeko"),
    ("langgraph", None, "langgraph>=0.2.0", "langgraph"),
    ("langchain_ollama", None, "langchain-ollama>=0.2.0", "langchain-ollama"),
    ("langchain_core", None, "langchain-core>=0.3.0", "langchain-core"),
    ("matplotlib", None, "matplotlib>=3.7.0", "matplotlib"),
    ("seaborn", None, "seaborn>=0.12.0", "seaborn"),
    ("scipy", None, "scipy>=1.10.0", "scipy"),
]

print("  패키지 점검 시작...")
failed = []

# numpy 버전 체크 (2.x면 다운그레이드)
if check_import("numpy") and not check_numpy_version():
    print("  numpy 2.x 감지 → 1.x 다운그레이드")
    pip_install("numpy>=1.24.0,<2.0.0")

for module, from_name, pip_pkg, desc in packages:
    if check_import(module, from_name):
        print(f"  [OK] {desc}")
    else:
        ok = pip_install(pip_pkg)
        # rdkit fallback: rdkit 실패 시 rdkit-pypi 시도
        if not ok and module == "rdkit":
            ok = pip_install("rdkit-pypi")
        if not ok:
            failed.append(desc)

if failed:
    print(f"\n  [WARNING] 설치 실패 패키지: {', '.join(failed)}")
    print("  파이프라인은 fallback으로 계속 진행합니다.")
else:
    print("\n  모든 패키지 설치 확인 완료!")
