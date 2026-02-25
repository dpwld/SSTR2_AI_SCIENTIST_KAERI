#!/usr/bin/env python3
"""수용체 PDBQT 포맷 검증 스크립트."""
import sys

if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

try:
    with open("SSTR2_receptor.pdbqt", "r") as f:
        lines = [l for l in f if l.startswith("ATOM")]

    print(f"  ATOM lines: {len(lines)}")
    if lines:
        sample = lines[0].rstrip()
        print(f"  Sample: {sample[:70]}...")
        print(f"  Line length: {len(sample)}")
        charge_part = sample[66:76].strip()
        type_part = sample[76:].strip()
        print(f"  Charge field: [{charge_part}]")
        print(f"  Type field: [{type_part}]")
    print("  검증 완료")
except FileNotFoundError:
    print("  [ERROR] SSTR2_receptor.pdbqt 파일을 찾을 수 없습니다.")
    sys.exit(1)
