#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  SSTR2 AI Scientist — 중앙 설정 파일 (config.py)
==========================================================================

  전체 파이프라인의 핵심 파라미터를 한 곳에서 관리한다.
  환경 변수(.env) 또는 OS 환경 변수를 통해 값을 오버라이드할 수 있어
  프로덕션(1억 개) 전환 시 코드 수정 없이 설정만 변경하면 된다.

  [사용 예시]
  테스트 모드 (기본):  python 00_peptide_gen.py
  프로덕션 모드:       NUM_SEQUENCES=100000000 python 00_peptide_gen.py
==========================================================================
"""

import os
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────
#  LLM 백엔드 설정 (후반전 Phase 2에서 사용)
# ─────────────────────────────────────────────────────────────────────────
OLLAMA_MODEL = os.getenv("OLLAMA_MODEL", "llama3.1:8b")

# ─────────────────────────────────────────────────────────────────────────
#  서열 생성 파라미터
# ─────────────────────────────────────────────────────────────────────────
# 테스트: 10,000 | 프로덕션: 100,000,000 (1억)
NUM_SEQUENCES = int(os.getenv("NUM_SEQUENCES", "10000"))

# 한 번에 메모리에 올리는 서열 수 (메모리 폭주 방지용 청크 크기)
CHUNK_SIZE = int(os.getenv("CHUNK_SIZE", "5000"))

# 병렬 워커 수 (CPU 코어 수에 맞춤)
NUM_WORKERS = int(os.getenv("NUM_WORKERS", str(os.cpu_count() or 4)))

# 재현성을 위한 랜덤 시드
RANDOM_SEED = int(os.getenv("RANDOM_SEED", "42"))

# ─────────────────────────────────────────────────────────────────────────
#  출력 디렉토리
# ─────────────────────────────────────────────────────────────────────────
OUTPUT_DIR = Path(os.getenv("OUTPUT_DIR", "output"))

# ─────────────────────────────────────────────────────────────────────────
#  대조군 Dotatate 기준 서열
# ─────────────────────────────────────────────────────────────────────────
# Dotatate (DOTA-Tyr3-octreotate, 상품명 Lutathera®) 펩타이드 부분:
#   실제 서열: D-Phe¹-Cys²-Tyr³-D-Trp⁴-Lys⁵-Thr⁶-Cys⁷-Thr⁸
#   D-아미노산을 포함하지만, 전하·pI·소수성 같은 스칼라 특성 계산에는
#   L-아미노산 근사(FCYWKTCT)를 사용한다.
DOTATATE_SEQ = "FCYWKTCT"

# ─────────────────────────────────────────────────────────────────────────
#  PK 필터 임계값 (레거시 — Module 01 5대 지표에서 내부적으로 참조)
# ─────────────────────────────────────────────────────────────────────────
# Net Charge 상한 (pH 7.4 기준). +1 초과 시 신장 독성 위험으로 탈락
MAX_NET_CHARGE = float(os.getenv("MAX_NET_CHARGE", "1.0"))

# pI(등전점) 허용 마진: Dotatate pI + 이 값까지 허용
PI_MARGIN = float(os.getenv("PI_MARGIN", "1.0"))

# Hydrophobicity(소수성) 허용 범위: Dotatate 값 ± 이 값
HYDRO_MARGIN = float(os.getenv("HYDRO_MARGIN", "0.5"))

# ─────────────────────────────────────────────────────────────────────────
#  Phase 1: 5대 안정성 지표 Composite Stability Score 가중치
# ─────────────────────────────────────────────────────────────────────────
# 각 지표의 정규화된 점수에 곱해지는 가중치.
# 합계 = 1.0. 값이 클수록 최종 순위에 더 큰 영향.
#
# [과학적 근거]
# (1) Protease Resistance (w=0.20): 절단 부위 수 → 혈중 반감기 직결
# (2) Aggregation Tendency (w=0.15): Boman index + KD → 제형 안정성
# (3) Kidney Toxicity (w=0.30): Net Charge → 신장 독성 (3대 지표 ③)
# (4) PPB (w=0.20): LogP + TPSA → 유리 약물 농도 → 표적 도달
# (5) Rigidity (w=0.15): Rotatable bonds → 결합 엔트로피 페널티
W_PROTEASE   = float(os.getenv("W_PROTEASE",   "0.20"))
W_AGGREGATION = float(os.getenv("W_AGGREGATION", "0.15"))
W_KIDNEY     = float(os.getenv("W_KIDNEY",     "0.30"))
W_PPB        = float(os.getenv("W_PPB",        "0.20"))
W_RIGIDITY   = float(os.getenv("W_RIGIDITY",   "0.15"))

# Module 01 출력 상위 후보 수 (Composite Score 기준 상위 N개)
FILTER_TOP_N = int(os.getenv("FILTER_TOP_N", "100"))

# ─────────────────────────────────────────────────────────────────────────
#  Phase 2: 3D 도킹 파이프라인 설정
# ─────────────────────────────────────────────────────────────────────────
# SSTR2 수용체 PDB/PDBQT 파일 경로 (로컬 고정)
RECEPTOR_PATH = Path(os.getenv("RECEPTOR_PATH", "SSTR2_receptor.pdbqt"))

# 도킹 엔진 실행파일 경로 (vina 또는 full path)
_default_vina = str(Path(__file__).resolve().parent / "vina_1.2.7_win.exe") \
    if (Path(__file__).resolve().parent / "vina_1.2.7_win.exe").exists() else "vina"
VINA_EXECUTABLE = os.getenv("VINA_EXECUTABLE", _default_vina)

# 도킹 검색 박스 중심 좌표 (Å) — SSTR2 orthosteric pocket
# PDB 7T10, Chain P (SST-14 리간드) CA 원자 중심으로부터 계산
DOCK_CENTER_X = float(os.getenv("DOCK_CENTER_X", "159.9"))
DOCK_CENTER_Y = float(os.getenv("DOCK_CENTER_Y", "149.6"))
DOCK_CENTER_Z = float(os.getenv("DOCK_CENTER_Z", "107.3"))

# 도킹 검색 박스 크기 (Å)
DOCK_SIZE_X = float(os.getenv("DOCK_SIZE_X", "30.0"))
DOCK_SIZE_Y = float(os.getenv("DOCK_SIZE_Y", "30.0"))
DOCK_SIZE_Z = float(os.getenv("DOCK_SIZE_Z", "30.0"))

# 도킹 exhaustiveness (Vina 정밀도, 높을수록 느리지만 정확)
DOCK_EXHAUSTIVENESS = int(os.getenv("DOCK_EXHAUSTIVENESS", "8"))

# ─────────────────────────────────────────────────────────────────────────
#  Phase 2: 대리 모델(Surrogate) 파라미터
# ─────────────────────────────────────────────────────────────────────────
# PeptiVerse 모델명 (HuggingFace). 펩타이드 특화 언어 모델 — 1순위 임베딩.
# ChatterjeeLab/PeptiVerse: 펩타이드 기능 예측에 최적화된 최신 PLM.
PEPTIVERSE_MODEL_NAME = os.getenv("PEPTIVERSE_MODEL_NAME", "ChatterjeeLab/PeptiVerse")

# ESM-2 모델명 (HuggingFace). PeptiVerse 실패 시 2순위 Fallback.
# 8M 파라미터 경량 모델이 기본값. GPU 가용 시 대형 모델로 교체 가능.
ESM2_MODEL_NAME = os.getenv("ESM2_MODEL_NAME", "facebook/esm2_t6_8M_UR50D")

# 정밀 도킹 시뮬레이션 대상 서열 수 (XGBoost 훈련용 레이블 생성)
TRAINING_SUBSET_SIZE = int(os.getenv("TRAINING_SUBSET_SIZE", "500"))

# 최종 AI 드림팀 원탁회의에 전달할 상위 후보 수
TOP_K_CANDIDATES = int(os.getenv("TOP_K_CANDIDATES", "100"))

# 도킹 시뮬레이션 재현성 시드
DOCKING_RANDOM_SEED = int(os.getenv("DOCKING_RANDOM_SEED", "42"))

# ─────────────────────────────────────────────────────────────────────────
#  Phase 2: AI 드림팀 원탁회의 파라미터
# ─────────────────────────────────────────────────────────────────────────
# 원탁회의에서 실제 평가할 상위 후보 수 (LLM 호출 비용 고려)
BOARDROOM_TOP_N = int(os.getenv("BOARDROOM_TOP_N", "5"))

# 기준 미달 시 설계 생성기의 최대 재설계 반복 횟수
MAX_REDESIGN_ITERATIONS = int(os.getenv("MAX_REDESIGN_ITERATIONS", "3"))

# Ollama 서버 주소
OLLAMA_BASE_URL = os.getenv("OLLAMA_BASE_URL", "http://localhost:11434")

# ─────────────────────────────────────────────────────────────────────────
#  물리 상수 (Physical Constants)
# ─────────────────────────────────────────────────────────────────────────
# R (Gas constant) in cal/(mol·K)
R_CAL = 1.987  # cal/(mol·K)
# 체온 T = 310.15 K (37 °C)
T_BODY = 310.15  # K
# RT at body temperature in kcal/mol
RT_BODY_KCAL = R_CAL * T_BODY / 1000.0  # ≈ 0.616 kcal/mol
