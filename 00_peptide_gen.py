#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  Module 00: 초거대 SSTR2 표적 펩타이드 서열 생성기  (00_peptide_gen.py)
==========================================================================

  ┌──────────────────────────────────────────────────────────────────┐
  │  프로젝트: SSTR2 AI Scientist — 차세대 방사성의약품(PRRT) 설계  │
  │  단계: 전반전 Phase 1 — 초고속 깔때기 파이프라인 (1/2)          │
  └──────────────────────────────────────────────────────────────────┘

  [프로젝트 배경]
  SSTR2(소마토스타틴 수용체 2형)는 신경내분비종양(NET) 세포 표면에
  정상 세포 대비 수십 배 과발현되는 G-단백질 결합 수용체(GPCR)이다.
  이 수용체를 정밀 타격하는 방사성의약품(PRRT: Peptide Receptor
  Radionuclide Therapy)은 종양에 방사선을 집중 전달하여 암세포를
  선택적으로 파괴하는 차세대 치료법이다.

  본 모듈은 SSTR2에 결합할 가능성이 있는 대규모 펩타이드 후보군을
  생물학적 제약 조건 하에 무작위로 생성한다.

  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  [3대 성공 지표와 본 모듈의 설계 관계]
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  ① 정밀타격 (Selectivity: 100배 이상 결합)
     → FWKT 약리활성 모티프를 모든 서열에 필수 고정.
       이 모티프는 소마토스타틴-14의 결합 활성 중심부에서 유래했으며,
       SSTR2 결합 포켓의 핵심 잔기들과 직접 상호작용하여
       나노몰(nM) 수준의 결합 친화도를 보장한다.
       비표적 수용체(SSTR1,3,4,5)는 결합 포켓 구조가 달라
       FWKT와의 상호작용이 약하므로 선택성이 확보된다.

  ② 생존력 (Retention: 72시간 이상 체류)
     → Cys-Cys 이황화 결합(-S-S-)으로 환형(Cyclic) 구조를 형성.
       - 고리 구조는 엔트로피 페널티를 줄여 결합 자유에너지(ΔG) 유리
       - 베타-턴(β-turn) 구조를 안정화하여 해리 속도(k_off) 최소화
       - 프로테아제 접근을 물리적으로 차단하여 혈중 반감기 연장
     → 가변 영역에서 K/R 제거로 트립신 계열 효소의 절단점 제거.
       트립신은 K/R의 C-말단을 절단하는데, 이를 없앰으로써
       혈중에서의 효소 분해를 원천 차단한다.

  ③ 안전성 (Safety: 신장 독성 20% 개선)
     → K(Lys)/R(Arg) 제거로 펩타이드의 순 양전하를 최소화.
       양전하 펩타이드는 신장 근위세뇨관 세포 표면의 메갈린(Megalin)
       수용체와 정전기적으로 결합하여 세포 내로 재흡수된다.
       이 재흡수가 방사성 물질의 신장 축적을 유발하여 독성을 일으킨다.
       K/R을 제거하면 이 메커니즘 자체를 억제할 수 있다.

  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  [서열 구조 (Scaffold Architecture)]
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  전체 서열은 아래 고정 뼈대(Scaffold)를 따른다:

    [N-tail] ─ C ─ [Inner-L] ─ F─W─K─T ─ [Inner-R] ─ C ─ [C-tail]
    ╰─선형─╯   ╰═══════ 이황화 결합 고리 (Cyclic Ring) ═══════╯   ╰─선형─╯
                          FWKT 모티프는 반드시 고리 내부

  - 총 길이: 8~14 아미노산 (Truncation 전략)
  - 고정 6잔기: C, F, W, K, T, C
  - 가변 2~8잔기: N-tail, Inner-L, Inner-R, C-tail에 무작위 분배

  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  [FWKT 모티프의 분자 수준 작용 메커니즘]
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  - F(Phe, 페닐알라닌):
    수용체 결합 포켓의 방향족 잔기들과 π-π 스태킹(Stacking) 상호작용.
    소수성 접촉으로 결합 엔탈피(ΔH)에 기여.

  - W(Trp, 트립토판):
    'Master Key(마스터키)' — SSTR2 막관통 나선(TM-helix) 깊숙한
    소수성 포켓에 인돌(Indole) 고리가 삽입된다.
    이 상호작용은 SSTR 아형 간 선택성의 핵심 결정 인자이다.
    (참고: PDB 7T10, 7T11 — SSTR2-소마토스타틴 복합체 구조)

  - K(Lys, 라이신):
    수용체 TM3(3번 막관통 나선)의 보존된 Asp(3.32)와
    정전기적 염교(Salt Bridge)를 형성. 이 상호작용은 모든 SSTR
    아형에서 리간드 인식의 필수 요소이며, 결합 에너지의 핵심 기여자.
    ※ FWKT 내부의 이 K만 예외적으로 유지하는 이유.

  - T(Thr, 트레오닌):
    수용체 세포외 루프 2(ECL2)의 백본 원자들과 수소결합 네트워크를
    형성하여 결합 복합체를 추가적으로 안정화한다.

  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  [실행 방법]
  python 00_peptide_gen.py                            # 기본(10,000개)
  NUM_SEQUENCES=100000000 python 00_peptide_gen.py    # 1억 개

  [출력]
  output/SSTR2_Raw_Library.csv  (컬럼: Sequence, Length)
==========================================================================
"""

import numpy as np
import pandas as pd
import time
import sys
from multiprocessing import Pool

# ── 프로젝트 중앙 설정 로드 ──
from config import (
    NUM_SEQUENCES, CHUNK_SIZE, NUM_WORKERS,
    RANDOM_SEED, OUTPUT_DIR
)


# ==========================================================================
#  생물학적 상수 정의 (Biological Constants)
# ==========================================================================

# ── 약리활성 핵심 모티프 (Pharmacophore Motif) ──
# 소마토스타틴-14(SST-14: AGCKNFFWKTFTSC)의 결합 활성 중심부에서 유래.
# 이 4잔기 서열은 SSTR2 결합 포켓과 직접 상호작용하는 최소 약리활성체로,
# 20년 이상의 구조-활성 관계(SAR) 연구에서 그 필수성이 검증되었다.
# (참고: Reubi JC, PNAS 2000; Weckbecker G, Nat Rev Drug Discov 2003)
PHARMACOPHORE = "FWKT"

# ── 이황화 결합 형성 잔기 ──
# 두 Cysteine의 티올(-SH) 기가 산화적 조건에서 -S-S- 공유결합을 형성.
# 이 이황화 결합은 펩타이드를 고리(Cyclic) 구조로 제약(Constrain)한다.
#
# [구조적 의미]
# - 고리 구조는 자유도(Degree of Freedom)를 줄여 결합 시 엔트로피 손실을
#   감소시킴 → 결합 자유에너지(ΔG) 유리 → 결합 친화도(Kd) 향상 (지표①)
# - 프로테아제의 기질 인식에 필요한 펩타이드 연장 구조를 방해하여
#   효소 분해를 물리적으로 차단 → 혈중 안정성 향상 (지표②)
# - β-턴 구조를 유도하여 FWKT 모티프를 수용체 결합에 최적인
#   공간 배향(Spatial Orientation)으로 고정
DISULFIDE_RESIDUE = "C"

# ── 가변 영역 아미노산 풀 (Variable Residue Pool) ──
# 표준 20종 아미노산에서 아래 3종을 제외한 17종:
#
#   ✗ K(Lys, 라이신) 제외 이유:
#     - 트립신/트립신 유사 프로테아제의 절단 부위 → 효소 저항성 감소 (지표②)
#     - 생리적 pH에서 양전하(+1) → 신장 메갈린 매개 재흡수 증가 (지표③)
#     ※ 단, FWKT 모티프 내부의 K는 수용체 결합에 필수이므로 예외적으로 유지
#
#   ✗ R(Arg, 아르기닌) 제외 이유:
#     - K와 동일하게 트립신 절단 부위이며 양전하(+1)를 가짐 (지표②③)
#     - 구아니디늄(Guanidinium) 기는 K의 아미노기보다 더 강한 양전하를
#       띠어 신장 독성 위험이 K보다 높음
#
#   ✗ C(Cys, 시스테인) 제외 이유:
#     - 추가 Cys가 존재하면 의도치 않은 이황화 결합(Off-target S-S)이
#       형성될 수 있어 설계된 고리 구조가 파괴될 위험
#     - Cys 위치는 이미 고정 뼈대(Scaffold)에 2개가 확정되어 있음
#
# 최종 풀 (17종):
#   A(Ala), D(Asp), E(Glu), F(Phe), G(Gly), H(His), I(Ile),
#   L(Leu), M(Met), N(Asn), P(Pro), Q(Gln), S(Ser), T(Thr),
#   V(Val), W(Trp), Y(Tyr)
VARIABLE_AA_POOL = list("ADEFGHILMNPQSTVWY")

# ── 서열 길이 범위 (Truncation 전략) ──
# 최소 8잔기: 고정 6잔기(C+FWKT+C) + 최소 가변 2잔기
#   → 최소한의 구조적 다양성 확보, 신장 여과에 유리한 작은 크기 (지표③)
# 최대 14잔기: 고정 6잔기 + 최대 가변 8잔기
#   → 수용체 세포외 루프(ECL)와의 추가 접촉 면적 확보 (지표①)
#   → 너무 길면 면역원성(Immunogenicity) 증가 및 합성 비용 상승 위험
MIN_LENGTH = 8
MAX_LENGTH = 14

# 고정 잔기 수: C(1) + F(1) + W(1) + K(1) + T(1) + C(1) = 6
FIXED_RESIDUE_COUNT = 6


# ==========================================================================
#  핵심 생성 함수 (Core Generation Functions)
# ==========================================================================

def _random_partition(rng, n, k):
    """
    n개의 가변 잔기를 k개의 영역에 무작위 균일 분배한다.

    [알고리즘] Stars-and-Bars 확장
    - n개의 아이템을 k개의 구간으로 나누기 위해 (k-1)개의 구분점을
      [0, n] 범위에서 무작위 선택 후 정렬한다.
    - 인접 구분점 간의 차이가 각 영역의 크기가 된다.
    - 이 방법은 (0, 0, 0, n) 같은 편중 분배도 허용하여
      극도로 비대칭적인 구조(예: 긴 N-tail + 짧은 고리)도 탐색한다.

    [생물학적 의미]
    4개 영역(N-tail, Inner-L, Inner-R, C-tail)에 가변 잔기를
    다양하게 분배함으로써 서로 다른 3D 형태(Conformation)를 가진
    후보군을 최대한 넓게 탐색한다 (Conformational Diversity).
    - N-tail/C-tail이 긴 서열: 수용체 표면 ECL과 추가 접촉 가능
    - Inner가 긴 서열: 고리 내부 유연성 증가, FWKT 미세 배향 변화

    Parameters
    ----------
    rng : np.random.RandomState
        재현 가능한 워커별 독립 난수 생성기
    n : int
        분배할 가변 잔기 총 수 (2~8)
    k : int
        분배 대상 영역 수 (4: N-tail, Inner-L, Inner-R, C-tail)

    Returns
    -------
    list[int]
        각 영역에 할당된 잔기 수 리스트 (합 = n)
    """
    # (k-1)개의 구분점을 [0, n] 범위에서 균일 선택
    dividers = np.sort(rng.randint(0, n + 1, size=k - 1))

    # 구분점 사이의 간격 = 각 영역의 크기
    parts = []
    prev = 0
    for d in dividers:
        parts.append(int(d - prev))
        prev = d
    parts.append(int(n - prev))
    return parts


def _generate_batch(args):
    """
    멀티프로세싱 워커 함수: 한 배치(Batch)의 펩타이드 서열을 생성한다.

    [서열 조립 과정 — 5단계]
    각 서열은 다음 순서로 조립된다:
      1) 전체 길이 결정 (8~14, 균일 분포)
      2) 가변 잔기 수 계산 = 전체 길이 - 6 (고정 잔기)
      3) 가변 잔기를 4개 영역에 무작위 분배
      4) 각 영역의 잔기를 제한된 AA 풀(17종)에서 무작위 선택
      5) 최종 조립: [N-tail]-C-[Inner-L]-FWKT-[Inner-R]-C-[C-tail]

    [성능 최적화 전략]
    - 각 워커가 독립적인 RandomState를 사용 → 락(Lock) 경합 zero
    - NumPy 배열 인덱싱으로 아미노산 선택 → Python 루프 최소화
    - 문자열 리스트 → join()으로 한 번에 결합 → 메모리 할당 최소화
    - 1억 개 생성 시: NUM_WORKERS=16이면 워커당 ~625만 개 처리

    Parameters
    ----------
    args : tuple (batch_size, seed)
        batch_size : int  — 이 워커가 생성할 서열 수
        seed : int        — 재현성을 위한 워커별 고유 시드

    Returns
    -------
    list[str]
        생성된 펩타이드 서열 문자열 리스트
    """
    batch_size, seed = args

    # ── 워커별 독립 난수 생성기 초기화 ──
    # 각 워커가 고유한 시드를 가지므로 병렬 실행 시에도
    # 결과의 재현성(Reproducibility)이 보장된다.
    rng = np.random.RandomState(seed)

    # ── 아미노산 풀을 NumPy 배열로 변환 (인덱싱 최적화) ──
    aa_pool = np.array(VARIABLE_AA_POOL)
    pool_size = len(aa_pool)  # 17종

    results = []
    for _ in range(batch_size):
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        #  Step 1: 서열 길이 결정 (8~14 균일 분포)
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # - 짧은 서열(8~10): 신장 사구체 여과에 유리(안전성↑),
        #   합성 용이, 하지만 수용체 접촉 면적이 제한적
        # - 긴 서열(12~14): 수용체 ECL과 추가 접촉으로 친화도↑,
        #   하지만 면역원성 위험 및 합성 비용 증가
        # → 균일 분포로 모든 길이를 편향 없이 탐색
        length = rng.randint(MIN_LENGTH, MAX_LENGTH + 1)

        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        #  Step 2: 가변 잔기 수 계산
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # 전체 길이에서 고정 6잔기(C+FWKT+C)를 빼면 가변 잔기 수
        # 범위: 2 (길이 8일 때) ~ 8 (길이 14일 때)
        n_variable = length - FIXED_RESIDUE_COUNT

        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        #  Step 3: 4개 영역에 가변 잔기 분배
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # 분배 대상 영역:
        #   [0] N-tail    — 이황화 고리 밖 N-말단 연장부
        #   [1] Inner-L   — 고리 내부 좌측 (C와 FWKT 사이)
        #   [2] Inner-R   — 고리 내부 우측 (FWKT와 C 사이)
        #   [3] C-tail    — 이황화 고리 밖 C-말단 연장부
        parts = _random_partition(rng, n_variable, 4)
        n_tail, inner_left, inner_right, c_tail = parts

        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        #  Step 4: 가변 영역 아미노산 무작위 선택
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # 17종 풀에서 균일 분포로 선택 (K, R, C 완전 배제)
        # NumPy 인덱싱을 사용하여 Python random.choice보다 ~5배 빠름
        total_var = n_variable
        if total_var > 0:
            var_indices = rng.randint(0, pool_size, size=total_var)
            var_aas = aa_pool[var_indices]
        else:
            var_aas = np.array([], dtype="<U1")

        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        #  Step 5: 서열 조립 (Assembly)
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # 구조: [N-tail] - C - [Inner-L] - F-W-K-T - [Inner-R] - C - [C-tail]
        #       ╰선형 꼬리╯  ╰═══════ 이황화 고리 ═══════╯  ╰선형 꼬리╯
        seq_parts = []
        idx = 0

        # (a) N-tail: 이황화 고리 바깥 N-말단 연장부
        #     수용체 세포외 루프(ECL1, ECL2)와의 추가 접촉 가능.
        #     또한 DOTA 킬레이터를 N-말단에 부착할 때 스페이서 역할.
        if n_tail > 0:
            seq_parts.append("".join(var_aas[idx:idx + n_tail]))
            idx += n_tail

        # (b) 첫 번째 Cysteine: 이황화 결합 시작점
        #     이 잔기부터 고리 구조가 시작된다.
        seq_parts.append(DISULFIDE_RESIDUE)

        # (c) Inner-Left: 고리 내부 좌측 스페이서
        #     FWKT 모티프의 3D 배향(Orientation)에 직접 영향.
        #     이 영역이 길수록 FWKT의 유연성이 증가하여
        #     유도 적합(Induced Fit) 결합에 유리할 수 있다.
        if inner_left > 0:
            seq_parts.append("".join(var_aas[idx:idx + inner_left]))
            idx += inner_left

        # (d) FWKT 약리활성 모티프: 수용체 결합의 핵심 엔진
        #     이 4잔기는 절대 변형 불가(Immutable).
        seq_parts.append(PHARMACOPHORE)

        # (e) Inner-Right: 고리 내부 우측 스페이서
        #     결합 후 수용체 내재화(Internalization) 효율에 영향.
        #     내재화 효율이 높을수록 방사성 물질이 종양 세포 내부에
        #     오래 머물러 치료 효과가 극대화된다 (지표②).
        if inner_right > 0:
            seq_parts.append("".join(var_aas[idx:idx + inner_right]))
            idx += inner_right

        # (f) 두 번째 Cysteine: 이황화 결합 종결점
        #     첫 번째 Cys와 -S-S- 공유결합을 형성하여 고리를 닫는다.
        seq_parts.append(DISULFIDE_RESIDUE)

        # (g) C-tail: 이황화 고리 바깥 C-말단 연장부
        #     C-말단 변형(예: Thr-ol 아미드화)을 위한 연장 영역.
        #     약물동태(PK) 프로필 미세 조정에 활용 가능.
        if c_tail > 0:
            seq_parts.append("".join(var_aas[idx:idx + c_tail]))

        results.append("".join(seq_parts))

    return results


# ==========================================================================
#  메인 실행부 (Main Execution)
# ==========================================================================

def main():
    """
    SSTR2 표적 펩타이드 라이브러리를 대규모로 생성한다.

    [실행 흐름 — 5단계]
    1. 출력 디렉토리 생성
    2. 생성 작업을 NUM_WORKERS개의 워커에 균등 분배
    3. 멀티프로세싱 풀(Pool)로 병렬 생성 수행
    4. 중복 제거 (Deduplication)
    5. 생성된 서열을 청크(Chunk) 단위로 CSV 저장 (메모리 효율)

    [스케일링 전략]
    - 10,000개 (테스트): 단일 머신, 수 초 내 완료
    - 100,000,000개 (프로덕션): 멀티코어 병렬 + 청크 I/O
      → 16코어 기준 약 10~20분 예상
    """
    print("=" * 70)
    print("  SSTR2 AI Scientist — Module 00: 초거대 서열 생성기")
    print("=" * 70)
    print(f"  목표 서열 수     : {NUM_SEQUENCES:>15,} 개")
    print(f"  병렬 워커 수     : {NUM_WORKERS:>15} 개")
    print(f"  서열 길이 범위   : {MIN_LENGTH}~{MAX_LENGTH} 아미노산")
    print(f"  가변 영역 AA 풀  : {''.join(VARIABLE_AA_POOL)} ({len(VARIABLE_AA_POOL)}종)")
    print(f"  고정 모티프      : C-...-FWKT-...-C (이황화 고리 + 약리활성체)")
    print(f"  랜덤 시드        : {RANDOM_SEED}")
    print("=" * 70)

    start_time = time.time()

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    #  1단계: 출력 디렉토리 생성
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "SSTR2_Raw_Library.csv"

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    #  2단계: 워커별 작업량 분배
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 총 서열 수를 워커 수로 나누되, 나머지는 앞쪽 워커들에 1개씩 추가.
    # 예: 10,000개 / 4워커 = [2500, 2500, 2500, 2500]
    #     10,003개 / 4워커 = [2501, 2501, 2501, 2500]
    base_batch = NUM_SEQUENCES // NUM_WORKERS
    remainder = NUM_SEQUENCES % NUM_WORKERS

    # 마스터 RNG에서 워커별 고유 시드를 파생
    # → 마스터 시드가 같으면 워커 시드도 동일 → 완전한 재현성
    master_rng = np.random.RandomState(RANDOM_SEED)
    tasks = []
    for i in range(NUM_WORKERS):
        batch_size = base_batch + (1 if i < remainder else 0)
        worker_seed = master_rng.randint(0, 2**31)
        tasks.append((batch_size, worker_seed))

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    #  3단계: 멀티프로세싱 병렬 생성
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print(f"\n>>> 서열 생성 시작 (워커 {NUM_WORKERS}개 병렬 실행)...")

    all_sequences = []

    if NUM_WORKERS > 1:
        # 멀티프로세싱 풀 사용
        # imap_unordered: 먼저 완료된 워커의 결과부터 수거 (최대 효율)
        with Pool(processes=NUM_WORKERS) as pool:
            for batch_result in pool.imap_unordered(_generate_batch, tasks):
                all_sequences.extend(batch_result)
                progress = len(all_sequences) / NUM_SEQUENCES * 100
                print(
                    f"\r  진행률: {progress:6.1f}% "
                    f"({len(all_sequences):>10,}/{NUM_SEQUENCES:,})",
                    end="", flush=True,
                )
    else:
        # 단일 워커 모드 (디버깅 시 스택 트레이스가 깔끔)
        for task in tasks:
            batch_result = _generate_batch(task)
            all_sequences.extend(batch_result)

    gen_elapsed = time.time() - start_time
    print(f"\n  [완료] 생성: {len(all_sequences):,} 개 ({gen_elapsed:.2f}초)")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    #  4단계: 중복 제거 (Deduplication)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 서열 공간이 매우 넓어(최소 17^2 ≈ 289 ~ 최대 17^8 ≈ 69억 조합)
    # 10,000개 수준에서 중복 확률은 극히 낮지만, 데이터 무결성을 위해
    # 중복을 제거한다.
    # ※ 1억 개 스케일에서는 메모리 이슈가 있을 수 있으며,
    #    이 경우 Bloom Filter 등으로 대체 가능.
    print(f"\n>>> 중복 제거 중...")
    unique_sequences = list(dict.fromkeys(all_sequences))  # 순서 보존 dedup
    n_duplicates = len(all_sequences) - len(unique_sequences)

    if n_duplicates > 0:
        print(f"  [경고] 중복 {n_duplicates:,}개 제거 "
              f"-> 최종 {len(unique_sequences):,}개 고유 서열")
    else:
        print(f"  [완료] 중복 없음: {len(unique_sequences):,}개 고유 서열")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    #  5단계: CSV 저장 (청크 기반, 메모리 효율)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 1억 개 규모에서는 전체 DataFrame을 한 번에 만들면 RAM이 터지므로
    # 청크(Chunk) 단위로 나눠서 append 모드로 저장한다.
    # 예: 1억 개 × 14자 평균 ≈ 1.4GB 텍스트. 청크 5,000개씩이면
    #     20,000회 I/O로 안전하게 기록 가능.
    print(f"\n>>> CSV 저장 중: {output_path}")

    # 헤더 먼저 작성
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("Sequence,Length\n")

    # 청크 단위 append
    for i in range(0, len(unique_sequences), CHUNK_SIZE):
        chunk = unique_sequences[i:i + CHUNK_SIZE]
        df_chunk = pd.DataFrame({
            "Sequence": chunk,
            "Length": [len(seq) for seq in chunk],
        })
        df_chunk.to_csv(
            output_path, mode="a", header=False,
            index=False, encoding="utf-8",
        )
        progress = min(i + CHUNK_SIZE, len(unique_sequences))
        print(
            f"\r  저장 진행: {progress:>10,}/{len(unique_sequences):,}",
            end="", flush=True,
        )

    elapsed = time.time() - start_time

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    #  최종 통계 보고
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print(f"\n\n{'=' * 70}")
    print(f"  >>> 생성 완료!")
    print(f"  출력 파일  : {output_path}")
    print(f"  총 서열 수 : {len(unique_sequences):,} 개")
    print(f"  소요 시간  : {elapsed:.2f} 초")
    print(f"  처리 속도  : {len(unique_sequences) / max(elapsed, 0.001):,.0f} 서열/초")
    print(f"{'=' * 70}")

    # ── 길이 분포 요약 ──
    # 8~14 균일 분포로 생성했으므로 각 길이의 비율이 ~14.3%에 가까워야 정상
    sample_seqs = unique_sequences[:min(10000, len(unique_sequences))]
    lengths = [len(seq) for seq in sample_seqs]
    print(f"\n  [길이 분포 (샘플 {len(sample_seqs):,}개 기준)]")
    for l in range(MIN_LENGTH, MAX_LENGTH + 1):
        count = lengths.count(l)
        pct = count / len(lengths) * 100
        bar = "#" * int(pct / 2)
        print(f"    {l:2d} AA: {count:5d} ({pct:5.1f}%) {bar}")

    # ── 샘플 서열 출력 ──
    print(f"\n  [샘플 서열 (상위 5개)]")
    for seq in unique_sequences[:5]:
        fwkt_pos = seq.find("FWKT")
        # C 위치 찾기 (이황화 결합 잔기)
        c_positions = [i for i, aa in enumerate(seq) if aa == "C"]
        print(
            f"    {seq:<20s} "
            f"(길이={len(seq):2d}, FWKT위치={fwkt_pos}, "
            f"S-S결합={c_positions})"
        )

    return output_path


if __name__ == "__main__":
    main()
