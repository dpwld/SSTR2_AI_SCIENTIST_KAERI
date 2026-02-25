#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  Module 01: 5대 안정성 지표 기반 Composite Stability Score 필터
  (01_fast_pk_filter.py)
==========================================================================

  [Stage]  Phase 1 — Step 2/2 (초고속 깔때기 파이프라인)
  [Input]  output/SSTR2_Raw_Library.csv   (Module 00 출력, ~10K 서열)
  [Output] output/SSTR2_Filtered_Top100.csv (Composite Score 상위 100개)

  [5대 안정성 지표 — Composite Stability Score]
  ──────────────────────────────────────────────────────────────────────
  모든 지표는 로컬에서 pip install만으로 즉시 계산 가능한 라이브러리만 사용.
  외부 서버/API/딥러닝 모델 호출 없음. No Hallucination.

  ① Protease Resistance (w=0.20)
     Python re 모듈로 Trypsin (K|R) 및 Chymotrypsin (F|W|Y|L) 절단부위 카운트.
     단, FWKT 모티프 내부의 K는 예외 처리 (보호된 약리활성 잔기).
     [근거] Werle M & Bernkop-Schnürch A, Amino Acids 2006; 30(4):351-367

  ② Aggregation Tendency (w=0.15)
     peptides 라이브러리의 boman() (단백질 결합 잠재력 지수) +
     hydrophobicity(scale='KyteDoolitle') (Kyte-Doolittle 소수성 스케일).
     Boman > 2.48: 높은 응집 경향 / KD hydro > 0.5: 과소수성 → 응집 위험.
     [근거] Boman HG, J Intern Med 2003; 254(3):197-215

  ③ Kidney Toxicity (w=0.30)  ★ 3대 성공 지표 ③ 직결
     peptides의 Henderson-Hasselbalch 기반 pH 7.4 순전하 계산.
     Dotatate (+0.707) 대비 개선 여부를 정량 평가.
     [근거] Vegt E et al., J Nucl Med 2010; 51(7):1049-1058

  ④ PPB — Plasma Protein Binding (w=0.20)
     RDKit으로 서열→SMILES→Mol 변환 후:
     - Crippen.MolLogP (LogP: 지용성 지표)
     - Descriptors.CalcTPSA (극성 표면적)
     높은 LogP + 낮은 TPSA → 비특이적 혈장 단백질 결합 ↑ → 유리 약물 ↓
     [근거] Hitchcock SA & Pennington LD, J Med Chem 2006; 49(26):7559-7583

  ⑤ Rigidity (w=0.15)
     RDKit CalcNumRotatableBonds로 회전가능 결합 수 계산.
     회전 결합 많을수록 결합 시 엔트로피 페널티 증가 → 친화도 손실.
     [근거] Veber DF et al., J Med Chem 2002; 45(12):2615-2623

  [Composite Stability Score]
  CSS = w1·S_protease + w2·S_aggregation + w3·S_kidney
      + w4·S_ppb + w5·S_rigidity
  각 S_i ∈ [0, 1] (MinMax 정규화, 높을수록 좋음)
  CSS가 높을수록 약물 안정성 우수 → 상위 100개 추출.
==========================================================================
"""

import re
import sys
import time
import numpy as np
import pandas as pd
from multiprocessing import Pool
from peptides import Peptide

# ── RDKit 조건부 임포트 ──
# Prompt_5 절대 원칙: pip install rdkit-pypi로 설치 가능.
# 없으면 PPB/Rigidity 지표를 수학적 근사로 Fallback.
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Crippen, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# ── 프로젝트 중앙 설정 ──
from config import (
    CHUNK_SIZE, NUM_WORKERS, OUTPUT_DIR, DOTATATE_SEQ,
    MAX_NET_CHARGE, PI_MARGIN, HYDRO_MARGIN,
    W_PROTEASE, W_AGGREGATION, W_KIDNEY, W_PPB, W_RIGIDITY,
    FILTER_TOP_N,
)

# Windows UTF-8 콘솔 호환
if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


# ==========================================================================
#  FWKT 모티프 내 K 위치 탐색 헬퍼
# ==========================================================================

def _find_fwkt_k_positions(sequence):
    """FWKT 모티프 내부의 K 위치(0-indexed)를 반환한다."""
    positions = set()
    start = 0
    while True:
        idx = sequence.find("FWKT", start)
        if idx < 0:
            break
        positions.add(idx + 2)  # K는 FWKT[2]
        start = idx + 1
    return positions


# ==========================================================================
#  지표 ① Protease Resistance
# ==========================================================================

def calc_protease_resistance(sequence):
    """
    Trypsin/Chymotrypsin 절단 부위 수를 카운트하여 단백질 분해 저항성을 평가.

    - Trypsin 절단: K, R 뒤 (C-terminal side) — 단, FWKT 내부 K 제외
    - Chymotrypsin 절단: F, W, Y, L 뒤
    - FWKT 모티프 내 F, W는 보호된 약리활성 잔기이므로 절단 대상에서 제외

    Returns: int (절단 부위 수, 낮을수록 좋음)
    """
    protected_k = _find_fwkt_k_positions(sequence)
    fwkt_start = sequence.find("FWKT")

    cleavage_count = 0
    for i, aa in enumerate(sequence):
        if i == len(sequence) - 1:
            continue  # C-말단 잔기는 절단 불가 (뒤에 잔기 없음)
        # Trypsin: K, R
        if aa in "KR":
            if aa == "K" and i in protected_k:
                continue  # FWKT 내 K는 보호
            cleavage_count += 1
        # Chymotrypsin: F, W, Y, L
        elif aa in "FWYL":
            # FWKT 모티프 내부의 F, W는 보호
            if fwkt_start >= 0 and fwkt_start <= i <= fwkt_start + 3:
                continue
            cleavage_count += 1

    return cleavage_count


# ==========================================================================
#  지표 ② Aggregation Tendency
# ==========================================================================

def calc_aggregation_tendency(sequence):
    """
    Boman Index와 Kyte-Doolittle 소수성으로 응집 경향을 수치화.

    - Boman Index: 단백질 표면 잔류 에너지 합 (kcal/mol).
      값이 높을수록 다른 단백질과의 비특이적 결합 경향 ↑
    - Kyte-Doolittle Hydrophobicity: 양수 = 소수성, 음수 = 친수성.
      과도한 소수성 → 수용액 내 응집 위험.

    Returns: float (응집 점수, 낮을수록 좋음)
    """
    pep = Peptide(sequence)
    boman = pep.boman()
    kd_hydro = pep.hydrophobicity("KyteDoolittle")

    # 선형 결합: Boman이 높거나 KD 소수성이 높으면 응집 위험 증가
    # Boman ∈ [-3, 5] 범위, KD ∈ [-2, 2] 범위
    agg_score = max(boman, 0) + max(kd_hydro, 0)
    return round(agg_score, 4)


# ==========================================================================
#  지표 ③ Kidney Toxicity
# ==========================================================================

def calc_kidney_toxicity(sequence):
    """
    pH 7.4에서의 Henderson-Hasselbalch 기반 순전하.
    Dotatate (+0.707)보다 낮을수록 신장 독성 개선.

    Returns: float (순전하, 낮을수록 좋음)
    """
    pep = Peptide(sequence)
    return round(pep.charge(pH=7.4), 4)


# ==========================================================================
#  지표 ④ PPB (Plasma Protein Binding)
# ==========================================================================

def _seq_to_mol(sequence):
    """
    아미노산 서열을 RDKit Mol 객체로 변환한다.
    MolFromSequence (peptide mode) → SMILES → Mol 경로를 시도한다.
    """
    if not RDKIT_AVAILABLE:
        return None
    try:
        # RDKit의 MolFromSequence는 펩타이드 HELM/단일문자 서열을 지원
        mol = Chem.MolFromSequence(sequence)
        if mol is not None:
            return mol
    except Exception:
        pass
    try:
        # Fallback: MolFromFASTA
        fasta = f">peptide\n{sequence}"
        mol = Chem.MolFromFASTA(fasta)
        if mol is not None:
            return mol
    except Exception:
        pass
    return None


def calc_ppb_score(sequence):
    """
    RDKit으로 LogP(지용성)와 TPSA(극성 표면적)를 계산하여
    혈장 단백질 비특이적 결합률을 추정.

    - 높은 LogP + 낮은 TPSA → PPB ↑ → 유리 약물 ↓ → 표적 도달 저하
    - PPB 추정 점수 = LogP - 0.01 * TPSA (경험적 스케일링)

    [RDKit 미설치 시 Fallback]
    아미노산 소수성 지수 합 / TPSA 근사값 (잔기별 극성 기여도 합산)

    Returns: float (PPB 추정 점수, 낮을수록 좋음)
    """
    mol = _seq_to_mol(sequence)
    if mol is not None:
        try:
            logp = Crippen.MolLogP(mol)
            tpsa = Descriptors.CalcTPSA(mol)
            ppb_score = logp - 0.01 * tpsa
            return round(ppb_score, 4)
        except Exception:
            pass

    # Fallback: 아미노산 기반 수학적 근사
    # LogP 근사: Kyte-Doolittle hydrophobicity 합
    # TPSA 근사: 극성 잔기(D, E, K, R, H, N, Q, S, T, Y, C) 기여
    pep = Peptide(sequence)
    kd_hydro = pep.hydrophobicity("KyteDoolittle")
    logp_approx = kd_hydro * len(sequence)

    polar_aa = set("DEKRHNQSTYC")
    n_polar = sum(1 for aa in sequence if aa in polar_aa)
    tpsa_approx = n_polar * 40.0 + 60.0  # 대략적 TPSA 추정 (backbone + sidechains)

    ppb_score = logp_approx - 0.01 * tpsa_approx
    return round(ppb_score, 4)


# ==========================================================================
#  지표 ⑤ Rigidity (구조 강직성)
# ==========================================================================

def calc_rigidity(sequence):
    """
    RDKit CalcNumRotatableBonds로 회전 가능 결합 수 계산.
    회전 결합이 적을수록 결합 시 엔트로피 손실 감소 → 친화도 유리.

    [RDKit 미설치 시 Fallback]
    Pro(P) 수 = 강직화 잔기, Gly(G) 수 = 유연화 잔기로 근사.
    서열 길이에 비례하는 기본 회전 결합 수에서 Pro 보정.

    Returns: int (회전 가능 결합 수, 낮을수록 좋음)
    """
    mol = _seq_to_mol(sequence)
    if mol is not None:
        try:
            rot_bonds = Descriptors.CalcNumRotatableBonds(mol)
            return rot_bonds
        except Exception:
            pass

    # Fallback: 수학적 근사
    # 펩타이드 backbone: (n-1) 개의 펩타이드 결합 중 phi/psi 회전
    # Pro는 phi를 고정 → 회전 결합 -1
    # Cys-Cys 디설파이드: 루프 형성 → 회전 결합 -2
    n = len(sequence)
    base_rotatable = max(0, (n - 1) * 2)  # phi + psi per residue
    pro_count = sequence.count("P")
    cys_pairs = sequence.count("C") // 2
    estimated = max(0, base_rotatable - pro_count - cys_pairs * 4)
    return estimated


# ==========================================================================
#  Dotatate Baseline 5대 지표 계산
# ==========================================================================

def compute_dotatate_baseline():
    """Dotatate의 5대 안정성 지표를 계산하여 비교 기준으로 반환한다."""
    seq = DOTATATE_SEQ
    print(f"  대조군 서열 (L-AA 근사): {seq}")
    print(f"  실제 서열: D-Phe-Cys-Tyr-D-Trp-Lys-Thr-Cys-Thr")

    pep = Peptide(seq)
    baseline = {
        "protease": calc_protease_resistance(seq),
        "aggregation": calc_aggregation_tendency(seq),
        "charge": pep.charge(pH=7.4),
        "pI": pep.isoelectric_point(),
        "hydrophobicity_eisenberg": pep.hydrophobicity("Eisenberg"),
        "ppb": calc_ppb_score(seq),
        "rigidity": calc_rigidity(seq),
    }

    print(f"  ├─ ① Protease sites     : {baseline['protease']} (절단 부위)")
    print(f"  ├─ ② Aggregation score   : {baseline['aggregation']:.4f}")
    print(f"  ├─ ③ Net Charge (pH 7.4) : {baseline['charge']:+.4f}")
    print(f"  ├─ ④ PPB score           : {baseline['ppb']:.4f}")
    print(f"  └─ ⑤ Rigidity (rot.bonds): {baseline['rigidity']}")

    return baseline


# ==========================================================================
#  배치 계산 워커 (Multiprocessing Worker)
# ==========================================================================

def _compute_batch(sequences):
    """
    서열 배치에 대해 5대 안정성 지표 + 기본 물리화학 특성을 계산한다.
    """
    results = []
    for seq in sequences:
        try:
            pep = Peptide(seq)
            charge = pep.charge(pH=7.4)
            pi = pep.isoelectric_point()
            hydro_eis = pep.hydrophobicity("Eisenberg")

            protease = calc_protease_resistance(seq)
            aggregation = calc_aggregation_tendency(seq)
            ppb = calc_ppb_score(seq)
            rigidity = calc_rigidity(seq)

            results.append((
                seq, len(seq),
                charge, pi, hydro_eis,
                protease, aggregation, ppb, rigidity,
            ))
        except Exception:
            continue

    return results


# ==========================================================================
#  Composite Stability Score 산출
# ==========================================================================

def compute_composite_score(df, baseline):
    """
    5대 지표를 MinMax 정규화 후 가중 합산하여 Composite Stability Score를 산출.
    CSS ∈ [0, 1], 높을수록 안정성 우수.

    [정규화 방향]
    - Protease, Aggregation, Kidney(charge), PPB, Rigidity 모두
      '낮을수록 좋은' 지표이므로, 정규화 후 (1 - normalized)로 반전.
    """
    scores = pd.DataFrame(index=df.index)

    # ① Protease: 절단 부위 적을수록 좋음
    col = df["ProteaseSites"].values.astype(float)
    mn, mx = col.min(), col.max()
    if mx > mn:
        scores["S_protease"] = 1.0 - (col - mn) / (mx - mn)
    else:
        scores["S_protease"] = 1.0

    # ② Aggregation: 응집 점수 낮을수록 좋음
    col = df["AggregationScore"].values.astype(float)
    mn, mx = col.min(), col.max()
    if mx > mn:
        scores["S_aggregation"] = 1.0 - (col - mn) / (mx - mn)
    else:
        scores["S_aggregation"] = 1.0

    # ③ Kidney: 순전하의 절대값이 낮을수록 좋음 (음전하 ~ 중성이 최적)
    # Dotatate 전하(+0.707) 대비 개선도 반영
    col = df["NetCharge"].values.astype(float)
    # 양전하 페널티: charge > 0이면 그대로, charge < 0이면 0으로 처리
    kidney_penalty = np.maximum(col, 0.0)
    mn, mx = kidney_penalty.min(), kidney_penalty.max()
    if mx > mn:
        scores["S_kidney"] = 1.0 - (kidney_penalty - mn) / (mx - mn)
    else:
        scores["S_kidney"] = 1.0

    # ④ PPB: 점수 낮을수록 좋음 (유리 약물 농도 ↑)
    col = df["PPBScore"].values.astype(float)
    mn, mx = col.min(), col.max()
    if mx > mn:
        scores["S_ppb"] = 1.0 - (col - mn) / (mx - mn)
    else:
        scores["S_ppb"] = 1.0

    # ⑤ Rigidity: 회전 결합 적을수록 좋음
    col = df["Rigidity"].values.astype(float)
    mn, mx = col.min(), col.max()
    if mx > mn:
        scores["S_rigidity"] = 1.0 - (col - mn) / (mx - mn)
    else:
        scores["S_rigidity"] = 1.0

    # Composite Stability Score (가중 합산)
    css = (W_PROTEASE   * scores["S_protease"]
         + W_AGGREGATION * scores["S_aggregation"]
         + W_KIDNEY     * scores["S_kidney"]
         + W_PPB        * scores["S_ppb"]
         + W_RIGIDITY   * scores["S_rigidity"])

    return css, scores


# ==========================================================================
#  Hard Cutoff 필터 (안전성 최소 기준)
# ==========================================================================

def apply_hard_cutoffs(df, baseline):
    """
    CSS 순위 산출 전, 절대적 안전 기준 미달 서열을 탈락시킨다.
    - Net Charge > +1.0: 신장 독성 고위험 → 즉시 탈락
    - pI > Dotatate + margin: 과도한 염기성 → 탈락
    """
    total = len(df)
    pi_upper = baseline["pI"] + PI_MARGIN

    mask_charge = df["NetCharge"] <= MAX_NET_CHARGE
    mask_pi = df["pI"] <= pi_upper

    mask_all = mask_charge & mask_pi
    n_pass = mask_all.sum()
    n_fail = total - n_pass

    print(f"  [Hard Cutoff] NetCharge <= {MAX_NET_CHARGE:+.1f}: "
          f"{mask_charge.sum():,} 통과")
    print(f"  [Hard Cutoff] pI <= {pi_upper:.2f}: "
          f"{(mask_charge & mask_pi).sum():,} 통과")
    print(f"  Hard Cutoff 결과: {n_pass:,} / {total:,} 통과 "
          f"({n_fail:,} 탈락)")

    return df[mask_all].copy()


# ==========================================================================
#  Main Execution
# ==========================================================================

def main():
    print("=" * 70)
    print("  SSTR2 AI Scientist — Module 01: 5대 안정성 지표 필터")
    print("=" * 70)
    start_time = time.time()

    # 입출력 경로
    input_path = OUTPUT_DIR / "SSTR2_Raw_Library.csv"
    output_path = OUTPUT_DIR / "SSTR2_Filtered_Top100.csv"

    if not input_path.exists():
        print(f"\n  [오류] 입력 파일 없음: {input_path}")
        print(f"  → 먼저 python 00_peptide_gen.py를 실행하세요.")
        sys.exit(1)

    # ── Step 1: Dotatate Baseline ──
    print(f"\n>>> Step 1: 대조군(Dotatate) Baseline 5대 지표 계산")
    print(f"  RDKit: {'v' + Chem.rdBase.rdkitVersion if RDKIT_AVAILABLE else '미설치 (Fallback 모드)'}")
    baseline = compute_dotatate_baseline()

    # ── Step 2: 원본 로드 ──
    print(f"\n>>> Step 2: 원본 라이브러리 로드")
    df_raw = pd.read_csv(input_path, encoding="utf-8")
    total = len(df_raw)
    print(f"  로드: {total:,}개 서열")

    # ── Step 3: 5대 지표 병렬 계산 ──
    print(f"\n>>> Step 3: 5대 안정성 지표 병렬 계산 (워커 {NUM_WORKERS}개)")
    all_sequences = df_raw["Sequence"].tolist()

    chunk_size = max(1, len(all_sequences) // NUM_WORKERS)
    seq_chunks = [all_sequences[i:i+chunk_size]
                  for i in range(0, len(all_sequences), chunk_size)]

    calc_start = time.time()
    all_results = []

    if NUM_WORKERS > 1 and len(seq_chunks) > 1:
        with Pool(processes=NUM_WORKERS) as pool:
            for batch in pool.imap(_compute_batch, seq_chunks):
                all_results.extend(batch)
                pct = len(all_results) / total * 100
                print(f"\r  계산 진행: {pct:6.1f}% "
                      f"({len(all_results):>8,}/{total:,})",
                      end="", flush=True)
    else:
        for chunk in seq_chunks:
            all_results.extend(_compute_batch(chunk))

    calc_elapsed = time.time() - calc_start
    print(f"\n  [완료] {len(all_results):,}개 계산 ({calc_elapsed:.2f}초)")

    # DataFrame 구성
    df = pd.DataFrame(all_results, columns=[
        "Sequence", "Length",
        "NetCharge", "pI", "Hydrophobicity",
        "ProteaseSites", "AggregationScore", "PPBScore", "Rigidity",
    ])

    # ── Step 4: Hard Cutoff ──
    print(f"\n>>> Step 4: Hard Cutoff 필터 적용")
    df = apply_hard_cutoffs(df, baseline)

    # ── Step 5: Composite Stability Score ──
    print(f"\n>>> Step 5: Composite Stability Score 산출")
    print(f"  가중치: Protease={W_PROTEASE}, Aggregation={W_AGGREGATION}, "
          f"Kidney={W_KIDNEY}, PPB={W_PPB}, Rigidity={W_RIGIDITY}")

    css, sub_scores = compute_composite_score(df, baseline)
    df["CompositeStabilityScore"] = css.round(4)
    df["S_Protease"] = sub_scores["S_protease"].round(4)
    df["S_Aggregation"] = sub_scores["S_aggregation"].round(4)
    df["S_Kidney"] = sub_scores["S_kidney"].round(4)
    df["S_PPB"] = sub_scores["S_ppb"].round(4)
    df["S_Rigidity"] = sub_scores["S_rigidity"].round(4)

    # ── Step 6: 상위 N개 추출 ──
    top_n = min(FILTER_TOP_N, len(df))
    print(f"\n>>> Step 6: Composite Score 상위 {top_n}개 추출")

    df = df.sort_values("CompositeStabilityScore", ascending=False)
    df_top = df.head(top_n).copy()
    df_top.insert(0, "Rank", range(1, len(df_top) + 1))

    # ── Step 7: 저장 ──
    print(f"\n>>> Step 7: 결과 저장")
    df_top.to_csv(output_path, index=False, encoding="utf-8")
    print(f"  출력: {output_path}")
    print(f"  저장: {len(df_top):,}개 서열")

    elapsed = time.time() - start_time

    # ── 최종 보고서 ──
    print(f"\n{'=' * 70}")
    print(f"  >>> 5대 안정성 지표 필터링 완료!")
    print(f"  입력: {total:>10,}개 → 출력: {len(df_top):>6,}개")
    print(f"  소요 시간: {elapsed:.2f}초")
    print(f"{'=' * 70}")

    # 상위 5개 출력
    print(f"\n  [Top 5 후보 (CSS 기준)]")
    print(f"  {'Rank':>4} | {'Sequence':<16} | {'CSS':>6} | "
          f"{'Charge':>7} | {'Protease':>8} | {'PPB':>6} | {'Rigid':>5}")
    print(f"  {'-'*4}-+-{'-'*16}-+-{'-'*6}-+-{'-'*7}-+-{'-'*8}-+-{'-'*6}-+-{'-'*5}")
    for _, row in df_top.head(5).iterrows():
        print(f"  {int(row['Rank']):>4} | {row['Sequence']:<16} | "
              f"{row['CompositeStabilityScore']:>6.4f} | "
              f"{row['NetCharge']:>+6.3f} | "
              f"{int(row['ProteaseSites']):>8} | "
              f"{row['PPBScore']:>6.2f} | "
              f"{int(row['Rigidity']):>5}")

    # Dotatate 대비 개선도
    if len(df_top) > 0:
        avg_charge = df_top["NetCharge"].mean()
        avg_css = df_top["CompositeStabilityScore"].mean()
        print(f"\n  [Dotatate 대비]")
        print(f"  통과군 평균 전하: {avg_charge:+.4f} "
              f"(Dotatate: {baseline['charge']:+.4f})")
        if baseline["charge"] > 0:
            improvement = (1 - avg_charge / baseline["charge"]) * 100
            print(f"  순전하 개선율: {improvement:+.1f}%")
        print(f"  통과군 평균 CSS: {avg_css:.4f}")

    print(f"\n  [다음 단계]")
    print(f"  → 02_train_surrogate.py: 3D 도킹 + AI 대리 모델 학습")

    return output_path


if __name__ == "__main__":
    main()
