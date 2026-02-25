#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  SSTR2 수용체 PDBQT 파일 자동 준비 스크립트  (prepare_receptor.py)
==========================================================================

  PDB 7T10 (SSTR2-소마토스타틴 복합체)에서 수용체 구조를 추출하고
  AutoDock Vina용 PDBQT 형식으로 변환한다.

  [사용법]
    python prepare_receptor.py

  [출력]
    SSTR2_receptor.pdbqt   (프로젝트 루트)

  [도킹 박스 좌표]
    PDB 7T10의 리간드(SST-14) 중심 좌표를 자동 계산하여
    config.py 도킹 박스 설정을 출력한다.
==========================================================================
"""

import os
import sys
import urllib.request
import gzip
import io
from pathlib import Path

# Windows UTF-8
if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

PROJECT_DIR = Path(__file__).resolve().parent
OUTPUT_PDBQT = PROJECT_DIR / "SSTR2_receptor.pdbqt"

# PDB 7T10: SSTR2 bound to SST-14 analog (Cryo-EM structure)
PDB_ID = "7T10"
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb.gz"
PDB_URL_FALLBACK = f"https://files.rcsb.org/download/{PDB_ID}.pdb"


# ─── Gasteiger 부분전하 근사 (원소별 기본값) ─────────────────────────────
# 정밀 도킹에는 ADFR Suite/MGLTools의 전하 계산을 권장하지만,
# 초기 스크리닝에는 원소별 기본값으로 충분하다.
DEFAULT_CHARGES = {
    "C": 0.050, "N": -0.350, "O": -0.400, "S": -0.100,
    "H": 0.150, "FE": 0.500, "ZN": 1.000, "CA": 1.000,
    "MG": 1.000, "MN": 1.000, "P": 0.250,
}

# AD4 원자 타입 매핑
AD4_TYPES = {
    "C": "C", "CA": "C", "CB": "C", "CG": "C", "CD": "C",
    "CE": "C", "CZ": "C", "CH": "C",
    "N": "N", "NA": "NA", "NE": "N", "NH": "N", "ND": "N", "NZ": "N",
    "O": "OA", "OD": "OA", "OE": "OA", "OG": "OA", "OH": "OA",
    "OXT": "OA",
    "S": "SA", "SD": "SA", "SG": "SA",
    "H": "HD",
    "FE": "FE", "ZN": "ZN", "CA": "C", "MG": "MG",
}


def download_pdb(pdb_id):
    """RCSB에서 PDB 파일을 다운로드한다."""
    print(f"  PDB {pdb_id} 다운로드 중...")

    # 1) gzip 압축 버전 시도
    try:
        print(f"    URL: {PDB_URL}")
        req = urllib.request.Request(PDB_URL, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            compressed = resp.read()
        with gzip.open(io.BytesIO(compressed), "rt") as f:
            content = f.read()
        print(f"    다운로드 완료 ({len(content):,} bytes)")
        return content
    except Exception as e:
        print(f"    gzip 실패: {e}")

    # 2) 비압축 버전 시도
    try:
        print(f"    Fallback URL: {PDB_URL_FALLBACK}")
        req = urllib.request.Request(PDB_URL_FALLBACK, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            content = resp.read().decode("utf-8")
        print(f"    다운로드 완료 ({len(content):,} bytes)")
        return content
    except Exception as e:
        print(f"    Fallback 실패: {e}")
        raise RuntimeError(
            f"PDB {pdb_id} 다운로드 실패.\n"
            f"수동 다운로드: https://www.rcsb.org/structure/{pdb_id}\n"
            f"다운로드한 파일을 {PROJECT_DIR / f'{pdb_id}.pdb'}에 배치하세요."
        )


def parse_pdb_receptor(pdb_content, receptor_chain="R"):
    """
    PDB에서 수용체 체인의 ATOM 레코드를 추출한다.

    PDB 7T10 체인 구성:
    - R: SSTR2 수용체 (target)
    - Gs 복합체, nanobody 등: 도킹에 불필요

    물/이온/리간드는 제거한다.
    """
    receptor_lines = []
    ligand_coords = []  # 리간드 좌표 (도킹 박스 중심 계산용)
    hetatm_ligand = []

    # 먼저 어떤 체인이 있는지 확인
    chains = set()
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") and len(line) >= 22:
            chains.add(line[21])

    print(f"  PDB 체인 목록: {sorted(chains)}")

    # 수용체 체인 자동 감지: 'R' 또는 첫 번째 체인
    if receptor_chain not in chains:
        # 7T10에서 SSTR2는 보통 chain R 또는 chain A
        for candidate in ["R", "A", "B"]:
            if candidate in chains:
                receptor_chain = candidate
                break
        else:
            receptor_chain = sorted(chains)[0] if chains else "A"

    print(f"  수용체 체인 선택: {receptor_chain}")

    # 비단백질 잔기 (제외 대상)
    EXCLUDE_RESNAMES = {
        "HOH", "WAT", "SOL",  # 물
        "NA", "CL", "K", "MG", "ZN", "CA", "FE", "MN",  # 이온
        "GOL", "EDO", "PEG", "DMS", "SO4", "PO4",  # 용매/첨가제
    }

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM"):
            chain = line[21] if len(line) > 21 else ""
            resname = line[17:20].strip()

            if chain == receptor_chain and resname not in EXCLUDE_RESNAMES:
                receptor_lines.append(line)

        elif line.startswith("HETATM"):
            resname = line[17:20].strip()
            chain = line[21] if len(line) > 21 else ""

            # 리간드 좌표 수집 (도킹 박스 중심 계산용)
            if chain != receptor_chain and resname not in EXCLUDE_RESNAMES:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ligand_coords.append((x, y, z))
                except (ValueError, IndexError):
                    pass

    print(f"  수용체 ATOM 수: {len(receptor_lines)}")
    print(f"  리간드 좌표 수: {len(ligand_coords)}")

    return receptor_lines, ligand_coords


def atom_to_pdbqt(atom_line):
    """
    PDB ATOM 라인을 PDBQT 형식으로 변환한다.

    PDBQT = PDB columns 1-66 + partial_charge + AD4_atom_type
    """
    # 원자 이름 추출
    atom_name = atom_line[12:16].strip()

    # 원소 기호 결정
    if len(atom_line) >= 78 and atom_line[76:78].strip():
        element = atom_line[76:78].strip().upper()
    else:
        # 원자 이름에서 추론
        element = ""
        for ch in atom_name:
            if ch.isalpha():
                element = ch.upper()
                break
        if not element:
            element = "C"

    # AD4 원자 타입
    ad4_type = AD4_TYPES.get(atom_name, element)
    if ad4_type not in ("C", "N", "NA", "OA", "SA", "HD", "FE", "ZN", "MG"):
        # 기본 매핑
        if element == "C":
            ad4_type = "C"
        elif element == "N":
            ad4_type = "N"
        elif element == "O":
            ad4_type = "OA"
        elif element == "S":
            ad4_type = "SA"
        elif element == "H":
            ad4_type = "HD"
        else:
            ad4_type = element

    # 부분 전하
    charge = DEFAULT_CHARGES.get(element, 0.000)

    # PDBQT 라인 조립 — Vina 1.2.x 호환 컬럼 정렬
    # Col 1-66: PDB, Col 67-76: partial charge (10자리), Col 77-78: AD4 type
    pdb_part = atom_line[:66].ljust(66)
    return f"{pdb_part}    {charge:+.4f} {ad4_type:<2s}"


def compute_docking_box(ligand_coords, padding=5.0):
    """리간드 좌표로부터 도킹 박스 중심과 크기를 계산한다."""
    if not ligand_coords:
        print("  [경고] 리간드 좌표 없음 — 수용체 중심 사용")
        return None

    import numpy as np
    coords = np.array(ligand_coords)
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)
    size = span + 2 * padding  # 패딩 추가

    # 최소 크기 보장
    size = np.maximum(size, 20.0)

    return {
        "center_x": round(float(center[0]), 1),
        "center_y": round(float(center[1]), 1),
        "center_z": round(float(center[2]), 1),
        "size_x": round(float(size[0]), 1),
        "size_y": round(float(size[1]), 1),
        "size_z": round(float(size[2]), 1),
    }


def main():
    print("=" * 60)
    print("  SSTR2 수용체 PDBQT 파일 자동 준비")
    print("=" * 60)

    # ── Step 1: PDB 다운로드 또는 로컬 파일 사용 ──
    local_pdb = PROJECT_DIR / f"{PDB_ID}.pdb"
    if local_pdb.exists():
        print(f"\n>>> Step 1: 로컬 PDB 파일 사용: {local_pdb}")
        with open(local_pdb, "r", encoding="utf-8", errors="ignore") as f:
            pdb_content = f.read()
    else:
        print(f"\n>>> Step 1: PDB {PDB_ID} 다운로드")
        pdb_content = download_pdb(PDB_ID)
        # 로컬 캐시 저장
        with open(local_pdb, "w", encoding="utf-8") as f:
            f.write(pdb_content)
        print(f"  로컬 캐시 저장: {local_pdb}")

    # ── Step 2: 수용체 추출 ──
    print(f"\n>>> Step 2: 수용체 체인 추출")
    receptor_lines, ligand_coords = parse_pdb_receptor(pdb_content)

    if not receptor_lines:
        print("  [오류] 수용체 ATOM이 추출되지 않았습니다.")
        print("  PDB 파일을 확인하세요.")
        sys.exit(1)

    # ── Step 3: PDBQT 변환 ──
    print(f"\n>>> Step 3: PDB → PDBQT 변환")
    pdbqt_lines = []
    for line in receptor_lines:
        pdbqt_lines.append(atom_to_pdbqt(line))

    with open(OUTPUT_PDBQT, "w", encoding="utf-8") as f:
        for line in pdbqt_lines:
            f.write(line + "\n")
        f.write("END\n")

    print(f"  출력: {OUTPUT_PDBQT}")
    print(f"  ATOM 수: {len(pdbqt_lines)}")
    file_size = OUTPUT_PDBQT.stat().st_size
    print(f"  파일 크기: {file_size / 1024:.1f} KB")

    # ── Step 4: 도킹 박스 계산 ──
    print(f"\n>>> Step 4: 도킹 박스 좌표 계산")
    box = compute_docking_box(ligand_coords)

    if box:
        print(f"  리간드 중심 좌표:")
        print(f"    center_x = {box['center_x']}")
        print(f"    center_y = {box['center_y']}")
        print(f"    center_z = {box['center_z']}")
        print(f"  박스 크기:")
        print(f"    size_x = {box['size_x']}")
        print(f"    size_y = {box['size_y']}")
        print(f"    size_z = {box['size_z']}")

        # config.py 업데이트 제안
        print(f"\n  [안내] config.py에 아래 값을 반영하면 도킹 정밀도가 향상됩니다:")
        print(f"    DOCK_CENTER_X = {box['center_x']}")
        print(f"    DOCK_CENTER_Y = {box['center_y']}")
        print(f"    DOCK_CENTER_Z = {box['center_z']}")
        print(f"    DOCK_SIZE_X   = {box['size_x']}")
        print(f"    DOCK_SIZE_Y   = {box['size_y']}")
        print(f"    DOCK_SIZE_Z   = {box['size_z']}")

        # config.py에 도킹 박스 좌표 자동 반영
        _update_config_docking_box(box)
        return box
    else:
        print("  도킹 박스를 수동으로 설정하세요.")
        return None

    print(f"\n{'=' * 60}")
    print(f"  준비 완료! 이제 파이프라인을 실행할 수 있습니다.")
    print(f"{'=' * 60}")


def _update_config_docking_box(box):
    """config.py의 도킹 박스 좌표를 자동 업데이트한다."""
    config_path = PROJECT_DIR / "config.py"
    if not config_path.exists():
        print("  [경고] config.py를 찾을 수 없어 수동 업데이트가 필요합니다.")
        return

    content = config_path.read_text(encoding="utf-8")

    import re
    replacements = {
        "DOCK_CENTER_X": str(box["center_x"]),
        "DOCK_CENTER_Y": str(box["center_y"]),
        "DOCK_CENTER_Z": str(box["center_z"]),
        "DOCK_SIZE_X": str(box["size_x"]),
        "DOCK_SIZE_Y": str(box["size_y"]),
        "DOCK_SIZE_Z": str(box["size_z"]),
    }

    for var, val in replacements.items():
        # 패턴: DOCK_CENTER_X = float(os.getenv("DOCK_CENTER_X", "0.0"))
        pattern = rf'({var}\s*=\s*float\(os\.getenv\("{var}",\s*")([^"]*)(\")'
        replacement = rf'\g<1>{val}\3'
        content, count = re.subn(pattern, replacement, content)
        if count > 0:
            print(f"  config.py 업데이트: {var} = {val}")

    config_path.write_text(content, encoding="utf-8")
    print(f"  config.py 저장 완료")


if __name__ == "__main__":
    box = main()
