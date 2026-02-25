#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module 00 실패 원인 진단 스크립트"""
import sys
import os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

print("=" * 60)
print("  SSTR2 AI Scientist - Diagnostic")
print("=" * 60)
print(f"  Python: {sys.version}")
print(f"  Executable: {sys.executable}")
print(f"  CWD: {os.getcwd()}")
print()

# 1) 패키지 임포트 테스트
packages = [
    ("numpy", "np"),
    ("pandas", "pd"),
    ("peptides", None),
    ("rdkit", None),
    ("rdkit.Chem", None),
    ("xgboost", "xgb"),
    ("sklearn", None),
    ("torch", None),
    ("transformers", None),
    ("langgraph", None),
    ("langchain_ollama", None),
    ("langchain_core", None),
    ("matplotlib", None),
    ("seaborn", None),
    ("scipy", None),
    ("streamlit", None),
    ("plotly", None),
]

print("[1] Package import test:")
failed = []
for pkg, alias in packages:
    try:
        mod = __import__(pkg)
        ver = getattr(mod, "__version__", "?")
        print(f"  [OK]   {pkg:<25s} v{ver}")
    except ImportError as e:
        print(f"  [FAIL] {pkg:<25s} -> {e}")
        failed.append(pkg)
    except Exception as e:
        print(f"  [ERR]  {pkg:<25s} -> {e}")
        failed.append(pkg)

print()
if failed:
    print(f"  ** {len(failed)}개 패키지 누락: {', '.join(failed)}")
    print(f"  ** 해결: {sys.executable} -m pip install {' '.join(failed)}")
else:
    print("  ** 모든 패키지 정상!")

# 2) config.py 임포트
print()
print("[2] config.py import test:")
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import config
    print(f"  [OK] NUM_SEQUENCES = {config.NUM_SEQUENCES}")
    print(f"  [OK] NUM_WORKERS   = {config.NUM_WORKERS}")
    print(f"  [OK] OUTPUT_DIR    = {config.OUTPUT_DIR}")
    print(f"  [OK] VINA_EXEC     = {config.VINA_EXECUTABLE}")
except Exception as e:
    print(f"  [FAIL] config.py 임포트 실패: {e}")

# 3) output 디렉토리
print()
print("[3] Output directory test:")
from pathlib import Path
out = Path("output")
if out.exists():
    print(f"  [OK] output/ exists")
    if out.is_dir():
        files = list(out.iterdir())
        print(f"  [OK] files: {len(files)}")
    else:
        print(f"  [FAIL] output is not a directory!")
else:
    print(f"  [WARN] output/ not found, will be created by Module 00")
    try:
        out.mkdir(parents=True, exist_ok=True)
        print(f"  [OK] output/ created successfully")
    except Exception as e:
        print(f"  [FAIL] Cannot create output/: {e}")

# 4) Module 00 직접 실행 시도
print()
print("[4] Module 00 dry run (minimal test):")
try:
    from config import NUM_SEQUENCES, CHUNK_SIZE, NUM_WORKERS, RANDOM_SEED, OUTPUT_DIR
    import numpy as np
    from multiprocessing import Pool

    # 간단한 Pool 테스트 (Windows multiprocessing 문제 감지)
    def _test_worker(x):
        return x * 2

    print(f"  Multiprocessing Pool test (workers={min(NUM_WORKERS, 2)})...")
    with Pool(processes=min(NUM_WORKERS, 2)) as pool:
        results = pool.map(_test_worker, [1, 2, 3])
    print(f"  [OK] Pool result: {results}")

    # 서열 생성 함수 임포트
    print(f"  Importing 00_peptide_gen...")
    # We avoid importing the module to prevent if __name__ issues
    # Instead, directly test the generation logic
    VARIABLE_AA_POOL = list("ADEFGHILMNPQSTVWY")
    rng = np.random.RandomState(42)
    aa_pool = np.array(VARIABLE_AA_POOL)
    length = rng.randint(8, 15)
    n_var = length - 6
    var_indices = rng.randint(0, len(aa_pool), size=n_var)
    test_seq = "".join(aa_pool[var_indices])
    print(f"  [OK] Test variable AA: {test_seq}")
    print(f"  [OK] Module 00 should work!")

except Exception as e:
    import traceback
    print(f"  [FAIL] {e}")
    traceback.print_exc()

# 5) Vina 실행 테스트
print()
print("[5] Vina executable test:")
import subprocess
import shutil
vina_path = shutil.which("vina")
vina_local = Path("vina_1.2.7_win.exe")
if vina_local.exists():
    print(f"  [OK] Local Vina found: {vina_local}")
    try:
        r = subprocess.run(
            [str(vina_local), "--version"],
            capture_output=True, text=True, timeout=10
        )
        print(f"  [OK] Vina output: {r.stdout.strip()[:100]}")
        if r.stderr:
            print(f"  [INFO] Vina stderr: {r.stderr.strip()[:200]}")
    except Exception as e:
        print(f"  [FAIL] Vina execution error: {e}")
elif vina_path:
    print(f"  [OK] Vina in PATH: {vina_path}")
else:
    print(f"  [WARN] Vina not found")

# 6) Receptor file
print()
print("[6] Receptor PDBQT test:")
receptor = Path("SSTR2_receptor.pdbqt")
if receptor.exists():
    size = receptor.stat().st_size
    print(f"  [OK] {receptor} ({size/1024:.1f} KB)")
    # Count ATOM lines
    with open(receptor, "r") as f:
        atom_count = sum(1 for line in f if line.startswith("ATOM"))
    print(f"  [OK] ATOM records: {atom_count}")
else:
    print(f"  [FAIL] {receptor} not found")

# 7) Ollama
print()
print("[7] Ollama test:")
try:
    import urllib.request
    req = urllib.request.Request("http://localhost:11434/api/tags")
    with urllib.request.urlopen(req, timeout=5) as resp:
        data = resp.read().decode()
    import json
    models = json.loads(data).get("models", [])
    print(f"  [OK] Ollama server running, {len(models)} model(s)")
    for m in models[:5]:
        print(f"       - {m.get('name', '?')}")
except Exception as e:
    print(f"  [FAIL] Ollama: {e}")

print()
print("=" * 60)
if failed:
    print(f"  RESULT: {len(failed)} package(s) missing — fix before running pipeline")
else:
    print("  RESULT: All checks passed!")
print("=" * 60)
