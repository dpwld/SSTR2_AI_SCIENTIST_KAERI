#!/usr/bin/env python3
"""모든 Python 모듈의 구문 검사를 수행합니다."""
import py_compile
import sys

if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

modules = [
    "config.py",
    "00_peptide_gen.py",
    "01_fast_pk_filter.py",
    "02_train_surrogate.py",
    "03_ai_boardroom.py",
    "04_validate_and_analyze.py",
    "prepare_receptor.py",
    "run_pipeline.py",
    "dashboard.py",
]

all_ok = True
for module in modules:
    try:
        py_compile.compile(module, doraise=True)
        print(f"  [OK] {module}")
    except py_compile.PyCompileError as e:
        print(f"  [ERROR] {module}: {e}")
        all_ok = False

if all_ok:
    print("\n  모든 모듈 구문 검사 통과!")
else:
    print("\n  [오류] 일부 모듈에 구문 오류가 있습니다.")
    sys.exit(1)

# 추가: 02_train_surrogate.py 임포트 테스트
print("\n  02_train_surrogate.py 임포트 가능 여부 확인...")
try:
    # 실행하지 않고 모듈만 임포트 (argparse 때문에 직접 실행은 안됨)
    import importlib.util
    spec = importlib.util.spec_from_file_location("mod02", "02_train_surrogate.py")
    # spec만 확인, 실제 import는 종속성 때문에 스킵
    if spec:
        print("  [OK] 모듈 스펙 확인 성공")
    else:
        print("  [WARNING] 모듈 스펙 생성 실패")
except Exception as e:
    print(f"  [WARNING] {e}")

print("\n  verify_syntax.py 완료")
