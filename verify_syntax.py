#!/usr/bin/env python3
"""프로젝트 루트의 모든 Python 파일을 실행 없이 구문 검사합니다."""
from pathlib import Path
import sys

if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

project_dir = Path(__file__).resolve().parent
modules = sorted(project_dir.glob("*.py"))

all_ok = True
for module in modules:
    try:
        compile(module.read_bytes(), str(module), "exec")
        print(f"  [OK] {module.name}")
    except (SyntaxError, OSError) as e:
        print(f"  [ERROR] {module.name}: {e}")
        all_ok = False

if all_ok:
    print(f"\n  모든 {len(modules)}개 모듈 구문 검사 통과!")
else:
    print("\n  [오류] 일부 모듈에 구문 오류가 있습니다.")
    sys.exit(1)
