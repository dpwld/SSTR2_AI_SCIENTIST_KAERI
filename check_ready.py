#!/usr/bin/env python3
"""실제 파이프라인 실행 전 환경 진단."""
import sys, os, subprocess
from pathlib import Path

if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

print("=" * 60)
print("  SSTR2 Pipeline - Environment Check")
print("=" * 60)

ROOT = Path(__file__).resolve().parent
errors = []
warnings = []

# 1. Python version
print(f"\n[1] Python: {sys.version}")
print(f"    Executable: {sys.executable}")
v = sys.version_info
if v.major == 3 and v.minor >= 11:
    print(f"    OK (3.{v.minor})")
else:
    warnings.append(f"Python 3.{v.minor} - 3.11+ recommended")

# 2. Core packages
print(f"\n[2] Package Check:")
packages = {
    "numpy":        "import numpy; print(numpy.__version__)",
    "pandas":       "import pandas; print(pandas.__version__)",
    "peptides":     "import peptides; print('OK')",
    "rdkit":        "from rdkit import Chem; from rdkit.Chem import AllChem, Crippen, Descriptors; print('OK')",
    "xgboost":      "import xgboost; print(xgboost.__version__)",
    "sklearn":      "import sklearn; print(sklearn.__version__)",
    "torch":        "import torch; print(torch.__version__)",
    "transformers": "import transformers; print(transformers.__version__)",
    "sentencepiece":"import sentencepiece; print('OK')",
    "meeko":        "from meeko import MoleculePreparation, PDBQTWriterLegacy; print('OK')",
    "langgraph":    "import langgraph; print('OK')",
    "langchain_ollama": "from langchain_ollama import ChatOllama; print('OK')",
    "matplotlib":   "import matplotlib; print(matplotlib.__version__)",
    "seaborn":      "import seaborn; print(seaborn.__version__)",
    "scipy":        "import scipy; print(scipy.__version__)",
}

for name, cmd in packages.items():
    try:
        r = subprocess.run([sys.executable, "-c", cmd],
                          capture_output=True, text=True, timeout=30)
        if r.returncode == 0:
            ver = r.stdout.strip()
            print(f"    [OK]  {name:20s} {ver}")
        else:
            err = r.stderr.strip().split('\n')[-1][:80]
            print(f"    [MISS] {name:20s} {err}")
            if name in ("numpy", "pandas", "peptides", "rdkit", "xgboost", "sklearn"):
                errors.append(f"{name} not installed (REQUIRED)")
            elif name == "meeko":
                warnings.append("meeko not installed (will use rigid docking fallback)")
            elif name == "sentencepiece":
                warnings.append("sentencepiece missing (PeptiVerse will fail, ESM-2 fallback)")
            else:
                warnings.append(f"{name} not installed (optional)")
    except Exception as e:
        print(f"    [ERR]  {name:20s} {e}")

# 3. NumPy version check
print(f"\n[3] NumPy Version:")
try:
    import numpy
    major = int(numpy.__version__.split(".")[0])
    if major >= 2:
        errors.append(f"NumPy {numpy.__version__} (2.x) - rdkit incompatible! Need <2.0")
        print(f"    [ERROR] NumPy {numpy.__version__} - need <2.0 for rdkit")
    else:
        print(f"    [OK] NumPy {numpy.__version__}")
except ImportError:
    errors.append("NumPy not installed")

# 4. Vina
print(f"\n[4] AutoDock Vina:")
vina_path = ROOT / "vina_1.2.7_win.exe"
if vina_path.exists():
    print(f"    [OK] Found: {vina_path.name}")
    try:
        r = subprocess.run([str(vina_path), "--version"],
                          capture_output=True, text=True, timeout=10)
        ver = (r.stdout.strip() or r.stderr.strip())[:100]
        print(f"    Version: {ver}")
    except Exception as e:
        print(f"    [WARNING] Cannot execute: {e}")
        warnings.append("Vina found but cannot execute")
else:
    print(f"    [MISS] vina_1.2.7_win.exe not in project root")
    errors.append("Vina not found - required for production docking")

# 5. Receptor
print(f"\n[5] Receptor PDBQT:")
receptor = ROOT / "SSTR2_receptor.pdbqt"
if receptor.exists():
    size = receptor.stat().st_size
    with open(receptor, "r") as f:
        atom_count = sum(1 for l in f if l.startswith("ATOM"))
    print(f"    [OK] {receptor.name} ({size/1024:.1f} KB, {atom_count} atoms)")
else:
    print(f"    [MISS] SSTR2_receptor.pdbqt not found")
    print(f"    Fix: python prepare_receptor.py")
    warnings.append("Receptor PDBQT missing - run prepare_receptor.py first")

# 6. Ollama
print(f"\n[6] Ollama Server:")
try:
    import urllib.request
    r = urllib.request.urlopen("http://localhost:11434/api/tags", timeout=3)
    data = r.read().decode()
    print(f"    [OK] Server running")
    if "llama3" in data.lower() or "llama" in data.lower():
        print(f"    Model available")
    else:
        warnings.append("Ollama running but llama3.1:8b may not be pulled")
except Exception as e:
    print(f"    [MISS] Not running ({e})")
    warnings.append("Ollama not running - Module 03 will use mock mode")

# 7. Output directory
print(f"\n[7] Output Directory:")
out = ROOT / "output"
if out.exists():
    files = list(out.glob("*"))
    print(f"    Exists with {len(files)} files")
else:
    print(f"    Will be created")

# Summary
print(f"\n{'=' * 60}")
if errors:
    print(f"  ERRORS ({len(errors)}) - Must fix before running:")
    for e in errors:
        print(f"    [X] {e}")
if warnings:
    print(f"\n  WARNINGS ({len(warnings)}) - Pipeline may still run:")
    for w in warnings:
        print(f"    [!] {w}")

if not errors:
    print(f"\n  >>> READY TO RUN! <<<")
    print(f"  Command: python run_pipeline.py")
else:
    print(f"\n  >>> FIX ERRORS FIRST <<<")
    print(f"  Missing packages: pip install <package>")

print(f"{'=' * 60}")
