#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  Module 03: 5인 AI 드림팀 원탁회의  (03_ai_boardroom.py)
==========================================================================

  [Stage]  Phase 2 - Step 2/2 (후반전 최종 모듈)
  [Input]  output/SSTR2_Top100_Boardroom_Ready.csv (Module 02 출력)
  [Output] output/Boardroom_Transcript.log          (토론 전문 기록)
           output/Final_Winning_Peptides_Report.md   (최종 보고서)

  [시스템 개요]
  5명의 AI 전문가가 LangGraph StateGraph 위에서 실시간 협업한다.
  각 에이전트는 Ollama(로컬 LLM)를 두뇌로 사용하며,
  MCP(공동 작업 노트)를 통해 의견을 교환한다.

  [에이전트 구성 및 역할]
  1. 증거 큐레이터 (Evidence Curator)
     - 논문/PDB 구조 데이터 기반 지식 제공
     - Trp-Lys 마스터키 결합 메커니즘 검증
  2. 설계 생성기 (Design Generator)
     - 반려 시 핵심 뼈대(Cys+FWKT) 유지하며 변이 서열 제안
  3. 구조 분석가 (Structure Analyst) -- 지표 1,2 검증
     - 100배 정밀타격: 도킹 에너지 압도적 차이 확인
     - 72시간 생존력: 베타-턴 유지, 낮은 k_off 심사
  4. PK 전문가 (PK Expert) -- 지표 3 검증
     - 전하/친수성 분석으로 신장 독성 20% 개선 검증
  5. 실험 계획가 (Experiment Planner)
     - DOTA N-term 부착 전략 및 In vivo 실험 계획서 작성

  [State Machine (순환 구조)]
    증거큐레이터 -> 구조분석가 -> [PASS] -> PK전문가 -> [PASS] -> 실험계획가
                                  [REJECT]              [REJECT]
                                     |                     |
                                     v                     v
                                  설계생성기 <-----------+
                                     |
                                     v
                                  구조분석가 (재평가)
    최대 재설계 반복 초과 시 -> 다음 후보로 이동
==========================================================================
"""

import sys
import time
import json
import datetime
import textwrap
import numpy as np
import pandas as pd
from pathlib import Path
from typing import TypedDict, Annotated
from operator import add

# -- 프로젝트 설정 --
from config import (
    OUTPUT_DIR, OLLAMA_MODEL, OLLAMA_BASE_URL,
    BOARDROOM_TOP_N, MAX_REDESIGN_ITERATIONS, DOTATATE_SEQ
)

# Windows UTF-8
if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# =========================================================================
#  의존성 조건부 임포트 (Graceful Degradation)
# =========================================================================
LANGGRAPH_OK = False
OLLAMA_OK = False

try:
    from langgraph.graph import StateGraph, START, END
    LANGGRAPH_OK = True
except ImportError:
    pass

try:
    from langchain_ollama import ChatOllama
    from langchain_core.messages import SystemMessage, HumanMessage
    OLLAMA_OK = True
except ImportError:
    try:
        from langchain_community.chat_models import ChatOllama
        from langchain_core.messages import SystemMessage, HumanMessage
        OLLAMA_OK = True
    except ImportError:
        pass


# =========================================================================
#  에이전트 시스템 프롬프트 (Agent System Prompts) -- 하드코딩
# =========================================================================
# 각 에이전트의 전문 분야와 판단 기준을 명시한 시스템 프롬프트.
# prompt_2.txt의 요구사항을 그대로 반영하여 하드코딩한다.

PROMPT_EVIDENCE_CURATOR = (
    "You are a senior structural biologist specializing in GPCR-peptide interactions. "
    "Your role is to provide literature-backed evidence for the team. "
    "Analyze the candidate sequence against the Trp-Lys master-key binding mechanism "
    "established in PDB 7T10 (SSTR2-somatostatin complex): "
    "(1) Trp indole ring insertion into the TM3-TM7 hydrophobic pocket, "
    "(2) Lys salt bridge with Asp3.32, "
    "(3) Disulfide ring forming a stable beta-turn. "
    "Be concise. Cite specific PDB residue interactions. Do NOT repeat yourself."
)

PROMPT_STRUCTURE_ANALYST = (
    "You are a biophysics expert evaluating peptide-receptor binding. "
    "Evaluate the candidate against TWO success criteria plus structural feasibility: "
    "[Criterion 1: 100x Selectivity] The Selectivity value (Kd_nontarget / Kd_SSTR2) "
    "must be >= 100x. Check the Selectivity field directly — do NOT invent your own formula. "
    "[Criterion 2: 72h Retention] Binding energy must be <= -7.0 kcal/mol, which ensures "
    "k_off < 3.86e-6 /s for 72-hour target residence. Check the SSTR2_BindingEnergy field directly. "
    "Also evaluate k_off using the Rigidity indicator from Module 01 — fewer rotatable bonds "
    "= lower conformational entropy penalty = slower dissociation = longer retention. "
    "[Steric Hindrance Check] The bulky DOTA chelator (~550 Da) is conjugated at the N-terminus. "
    "Assess whether the N-terminal residue has sufficient spatial clearance — if a large aromatic "
    "residue (W, F, Y) is at position 1 AND inside the Cys-Cys ring, flag steric risk. "
    "Cys at position 1 is fine (DOTA conjugates via the free amine, not the thiol). "
    "If BOTH criteria are met and no fatal steric issue, respond with VERDICT: PASS. "
    "If EITHER fails, respond with VERDICT: REJECT and state which criterion failed and why. "
    "Be precise and use the actual numerical values provided."
)

PROMPT_PK_EXPERT = (
    "You are a pharmacokinetics and toxicology expert. "
    "Evaluate [Criterion 3: 20% Kidney Safety Improvement] using the 5-indicator "
    "Composite Stability Score from Module 01, NOT just Net Charge alone. "
    "The five indicators are: "
    "(1) Protease Resistance — fewer cleavage sites = longer plasma half-life. "
    "(2) Aggregation Tendency — lower score = better formulation stability. "
    "(3) Kidney Toxicity (Net Charge at pH 7.4) — Dotatate reference: +0.707. "
    "(4) PPB (Plasma Protein Binding via LogP/TPSA) — moderate PPB is desirable "
    "(too high = sequestered, too low = rapid clearance). "
    "(5) Rigidity (rotatable bonds) — fewer = less entropic penalty, slower k_off. "
    "EVALUATION RULES: "
    "- CompositeStabilityScore >= 0.5: strong PK profile. "
    "- CompositeStabilityScore 0.3–0.5: acceptable with caveats. "
    "- CompositeStabilityScore < 0.3: REJECT (poor overall stability). "
    "- Additionally, if NetCharge > +0.8: REJECT regardless of composite score "
    "(excessive renal accumulation risk via megalin-mediated reuptake). "
    "- If NetCharge <= +0.707 AND CompositeStabilityScore >= 0.5: VERDICT: PASS (excellent). "
    "Assess each indicator individually, then give an integrated VERDICT: PASS or VERDICT: REJECT. "
    "Use the actual numerical values provided."
)

PROMPT_DESIGN_GENERATOR = (
    "You are a senior peptide medicinal chemist redesigning a rejected SSTR2-targeting peptide. "
    "Provide a STRUCTURED response with ALL of the following sections:\n\n"
    "[1. REJECTION DIAGNOSIS]\n"
    "Identify which criterion failed (Selectivity <100x / Binding >-7.0 / Charge >+0.8 / CSS <0.3) "
    "and explain the biophysical root cause in 1-2 sentences.\n\n"
    "[2. DESIGN STRATEGY]\n"
    "Explain your mutation strategy — which residue positions you will change and WHY, citing "
    "structure-activity relationships. For example:\n"
    "- Selectivity issue → add aromatic residues (F,W,Y) flanking FWKT to enhance hydrophobic packing in TM3-TM7 pocket\n"
    "- Weak binding → optimize ring size or add H-bond donors (N,Q,S,T) near the Cys bridge\n"
    "- High charge → introduce acidic residues (D,E) or replace polar residues with neutral ones\n"
    "- Low CSS → address the weakest sub-indicator (protease sites, aggregation, rigidity)\n\n"
    "[3. MUTATION DETAILS]\n"
    "List each mutation as: Position X: OldAA → NewAA (rationale)\n\n"
    "[4. NEW SEQUENCE]\n"
    "NEW_SEQUENCE: <your sequence in uppercase>\n\n"
    "CONSTRAINTS (must follow ALL): "
    "1. Keep the Cys-Cys disulfide ring intact. "
    "2. Keep the FWKT pharmacophore intact. "
    "3. Only mutate residues OUTSIDE Cys ring and FWKT. "
    "4. Do NOT use K, R, or C in variable positions. "
    "5. Length: 8-14 amino acids."
)

PROMPT_EXPERIMENT_PLANNER = (
    "You are a radiochemistry and preclinical study planner. "
    "Write a concise experimental protocol for the approved peptide candidate. Include: "
    "(1) DOTA-NHS coupling: conditions (pH, temperature, time, purification). "
    "(2) Lu-177 radiolabeling: protocol, radiochemical purity target. "
    "(3) In vitro binding: AR42J cell Kd measurement, selectivity panel (CHO-SSTR1/3/4/5). "
    "(4) In vivo biodistribution: AR42J xenograft nude mouse, time points (1h, 4h, 24h, 48h, 72h), "
    "organs measured (%ID/g), tumor-to-kidney ratio. "
    "(5) Success criteria: tumor uptake >10 %ID/g at 4h, <50% decrease at 72h, "
    "tumor/kidney ratio >1.2x vs Dotatate. "
    "Be specific with numbers and conditions. Do NOT pad with filler text."
)

# =========================================================================
#  가변 영역 아미노산 풀 (변이 생성용)
# =========================================================================
VARIABLE_AA_POOL = list("ADEFGHILMNPQSTVWY")


# =========================================================================
#  LLM 백엔드 (Ollama + Mock Fallback)
# =========================================================================

def _check_ollama_live():
    """Ollama 서버 접속 가능 여부를 확인한다."""
    try:
        import urllib.request
        url = f"{OLLAMA_BASE_URL}/api/tags"
        req = urllib.request.Request(url, method="GET")
        resp = urllib.request.urlopen(req, timeout=3)
        return resp.status == 200
    except Exception:
        return False


def _create_llm():
    """Ollama LLM 인스턴스를 생성한다."""
    if OLLAMA_OK and _check_ollama_live():
        print(f"  [LLM] Ollama 연결 성공 (모델: {OLLAMA_MODEL})")
        return ChatOllama(
            model=OLLAMA_MODEL,
            base_url=OLLAMA_BASE_URL,
            temperature=0.7,
            timeout=120,
        )
    else:
        print(f"  [LLM] Ollama 미접속 -> Mock 모드로 실행")
        return None


def _call_llm(llm, system_prompt, user_prompt, timeout_sec=120):
    """LLM을 호출하여 응답을 받는다. 실패/타임아웃 시 None을 반환하여 Mock으로 전환."""
    if llm is not None:
        try:
            import threading
            result_container = [None]
            error_container = [None]

            def _invoke():
                try:
                    msgs = [SystemMessage(content=system_prompt), HumanMessage(content=user_prompt)]
                    result_container[0] = llm.invoke(msgs).content
                except Exception as e:
                    error_container[0] = e

            t = threading.Thread(target=_invoke, daemon=True)
            t.start()
            t.join(timeout=timeout_sec)

            if t.is_alive():
                print(f"    [LLM 타임아웃: {timeout_sec}초 초과] -> Mock 응답 사용",
                      flush=True)
                return None
            if error_container[0]:
                print(f"    [LLM 오류: {error_container[0]}] -> Mock 응답 사용",
                      flush=True)
                return None
            return result_container[0]
        except Exception as e:
            print(f"    [LLM 오류: {e}] -> Mock 응답 사용", flush=True)
            return None
    return None  # Mock 모드 시그널


# =========================================================================
#  Mock 응답 생성기 (Ollama 미접속 시)
# =========================================================================

def _mock_evidence_curator(seq, data):
    c_pos = [i for i, a in enumerate(seq) if a == "C"]
    ring_info = f"Cys({c_pos[0]})-Cys({c_pos[1]})" if len(c_pos) >= 2 else "N/A"
    ring_size = (c_pos[1] - c_pos[0] - 1) if len(c_pos) >= 2 else 0
    return (
        f"[Evidence Curator]\n"
        f"Sequence {seq} contains the FWKT pharmacophore (Phe-Trp-Lys-Thr), the master-key "
        f"motif for SSTR2 binding. Per PDB 7T10 (SSTR2-somatostatin complex):\n"
        f"- Trp indole ring inserts into the TM3-TM7 hydrophobic pocket (Phe6.52/Phe7.35 stacking).\n"
        f"- Lys forms a critical salt bridge with Asp3.32 on TM3.\n"
        f"- Disulfide bond {ring_info} (ring size: {ring_size} residues) constrains the peptide "
        f"into a beta-turn, mimicking the bioactive conformation of octreotide.\n"
        f"This structural arrangement is consistent with Weckbecker et al. (Nat Rev Drug Discov, 2003)."
    )


def _mock_structure_analyst(seq, data):
    energy = data.get("SSTR2_BindingEnergy", -7.0)
    sel = data.get("Selectivity", 50.0)
    if energy <= -7.0 and sel >= 100.0:
        verdict = "PASS"
        reason = (
            f"Criterion 1 (100x Selectivity): Selectivity = {sel:.1f}x (>= 100x threshold). PASS.\n"
            f"Criterion 2 (72h Retention): Binding energy = {energy:.3f} kcal/mol (<= -7.0 threshold). "
            f"This corresponds to k_off < 3.86e-6 /s, sufficient for 72h target residence. PASS.\n"
            f"The disulfide-constrained beta-turn maintains rigidity, further reducing k_off."
        )
    elif energy <= -7.0 and sel < 100.0:
        verdict = "REJECT"
        reason = (
            f"Criterion 1 (100x Selectivity): Selectivity = {sel:.1f}x (< 100x threshold). FAIL.\n"
            f"Criterion 2 (72h Retention): Binding energy = {energy:.3f} kcal/mol. PASS.\n"
            f"Rejection reason: Insufficient selectivity over non-target SSTR subtypes."
        )
    elif energy > -7.0 and sel >= 100.0:
        verdict = "REJECT"
        reason = (
            f"Criterion 1 (100x Selectivity): Selectivity = {sel:.1f}x. PASS.\n"
            f"Criterion 2 (72h Retention): Binding energy = {energy:.3f} kcal/mol (> -7.0 threshold). FAIL.\n"
            f"Rejection reason: Binding too weak for 72h residence time."
        )
    else:
        verdict = "REJECT"
        reason = (
            f"Criterion 1 (100x Selectivity): Selectivity = {sel:.1f}x (< 100x). FAIL.\n"
            f"Criterion 2 (72h Retention): Binding energy = {energy:.3f} kcal/mol (> -7.0). FAIL.\n"
            f"Rejection reason: Both selectivity and retention criteria unmet."
        )
    return f"[Structure Analyst]\n{reason}\nVERDICT: {verdict}"


def _mock_pk_expert(seq, data):
    charge = data.get("NetCharge", 0.0)
    css = data.get("CompositeStabilityScore", 0.5)
    protease = data.get("ProteaseSites", "N/A")
    agg = data.get("AggregationScore", "N/A")
    ppb = data.get("PPBScore", "N/A")
    rigidity = data.get("Rigidity", "N/A")

    # Evaluate each indicator
    indicator_lines = (
        f"  1. Protease Cleavage Sites: {protease} — "
        f"{'low risk' if isinstance(protease, (int, float)) and protease <= 2 else 'moderate/high risk'}.\n"
        f"  2. Aggregation Score: {agg} — "
        f"{'acceptable' if isinstance(agg, (int, float)) and agg < 2.0 else 'elevated'}.\n"
        f"  3. Net Charge (pH 7.4): {charge:+.4f} vs Dotatate +0.707 — "
        f"{'improved' if charge <= 0.707 else 'comparable' if charge <= 0.8 else 'elevated risk'}.\n"
        f"  4. PPB Score: {ppb} — "
        f"{'moderate (favorable)' if isinstance(ppb, (int, float)) and 0.3 <= ppb <= 0.7 else 'review needed'}.\n"
        f"  5. Rigidity: {rigidity} rotatable bonds — "
        f"{'rigid (favorable)' if isinstance(rigidity, (int, float)) and rigidity <= 10 else 'flexible'}."
    )

    if charge > 0.8:
        verdict = "REJECT"
        reason = (
            f"5-Indicator Stability Assessment:\n{indicator_lines}\n\n"
            f"Composite Stability Score: {css:.3f}\n"
            f"CRITICAL: Net Charge {charge:+.4f} exceeds +0.8 safety threshold.\n"
            f"Excessive positive charge drives megalin-mediated renal reabsorption.\n"
            f"Recommendation: introduce acidic residues (D, E) to reduce net charge."
        )
    elif css < 0.3:
        verdict = "REJECT"
        reason = (
            f"5-Indicator Stability Assessment:\n{indicator_lines}\n\n"
            f"Composite Stability Score: {css:.3f} (below 0.3 threshold).\n"
            f"Poor overall PK profile — multiple indicators below acceptable range."
        )
    elif charge <= 0.707 and css >= 0.5:
        verdict = "PASS"
        reason = (
            f"5-Indicator Stability Assessment:\n{indicator_lines}\n\n"
            f"Composite Stability Score: {css:.3f} (excellent, >= 0.5).\n"
            f"Net Charge {charge:+.4f} <= Dotatate (+0.707): kidney safety IMPROVED.\n"
            f"Strong overall PK profile across all five indicators."
        )
    else:
        verdict = "PASS"
        reason = (
            f"5-Indicator Stability Assessment:\n{indicator_lines}\n\n"
            f"Composite Stability Score: {css:.3f} (acceptable).\n"
            f"Net Charge {charge:+.4f} is within tolerable range.\n"
            f"Lysine co-infusion protocol recommended as standard PRRT practice."
        )
    return f"[PK Expert]\n{reason}\nVERDICT: {verdict}"


def _mock_design_generator(seq, rejection_reason):
    """
    반려 사유를 분석하여 과학적 근거에 기반한 변이를 제안한다.
    Cys ring + FWKT 뼈대는 절대 보존.
    """
    rng = np.random.RandomState(hash(seq + rejection_reason) % (2**31))
    seq_list = list(seq)

    # 보존 위치 식별
    c_pos = [i for i, a in enumerate(seq) if a == "C"]
    immutable = set(c_pos)
    fwkt_start = seq.find("FWKT")
    if fwkt_start >= 0:
        immutable.update(range(fwkt_start, fwkt_start + 4))
    mutable = [i for i in range(len(seq)) if i not in immutable]

    reason_upper = rejection_reason.upper()

    # ── 반려 사유 진단 ──
    diagnosis = ""
    strategy = ""
    mutations = []

    if "SELECTIVITY" in reason_upper and ("FAIL" in reason_upper or "< 100" in reason_upper):
        # 선택성 부족 → 소수성/방향족 잔기 강화 (SSTR2 TM3-TM7 포켓 패킹)
        diagnosis = (
            "Selectivity below 100x threshold. The peptide lacks sufficient hydrophobic "
            "complementarity with the SSTR2 TM3-TM7 pocket, allowing cross-reactivity "
            "with SSTR3/5 subtypes that have shallower binding grooves."
        )
        strategy = (
            "Introduce aromatic (F, Y, W) or bulky hydrophobic (I, L) residues flanking "
            "the FWKT motif to enhance van der Waals contacts with Phe6.52 and Phe7.35 "
            "in the SSTR2 pocket. These bulky sidechains create steric clashes in the "
            "narrower pockets of SSTR1/3/4/5, improving selectivity."
        )
        aromatic_pool = list("FWYIL")
        for pos in mutable:
            if seq_list[pos] not in aromatic_pool and rng.random() < 0.5:
                old_aa = seq_list[pos]
                new_aa = aromatic_pool[rng.randint(len(aromatic_pool))]
                seq_list[pos] = new_aa
                mutations.append(
                    f"Position {pos+1}: {old_aa} → {new_aa} "
                    f"(aromatic/bulky — enhances TM pocket hydrophobic packing)"
                )
                if len(mutations) >= 2:
                    break

    elif "BINDING" in reason_upper or "RETENTION" in reason_upper or "> -7.0" in reason_upper:
        # 결합 에너지 부족 → H-bond donor/acceptor 추가, 구조 강성 증가
        diagnosis = (
            "Binding energy weaker than -7.0 kcal/mol threshold. Insufficient hydrogen bonding "
            "network and/or conformational flexibility causing entropic penalty upon binding."
        )
        strategy = (
            "Replace non-polar residues at variable positions with H-bond capable residues "
            "(N, Q, S, T) to form additional contacts with SSTR2 ECL2 backbone carbonyls. "
            "Proline (P) can be introduced to reduce backbone flexibility, lowering the "
            "conformational entropy penalty (ΔS_conf) upon receptor binding."
        )
        hbond_pool = list("NQSTP")
        for pos in mutable:
            if seq_list[pos] not in hbond_pool and rng.random() < 0.5:
                old_aa = seq_list[pos]
                new_aa = hbond_pool[rng.randint(len(hbond_pool))]
                seq_list[pos] = new_aa
                rationale = ("Pro — restricts phi angle, reduces backbone flexibility"
                             if new_aa == "P" else
                             "H-bond donor/acceptor — additional polar contacts with ECL2")
                mutations.append(f"Position {pos+1}: {old_aa} → {new_aa} ({rationale})")
                if len(mutations) >= 2:
                    break

    elif "CHARGE" in reason_upper or "MEGALIN" in reason_upper or "+0.8" in reason_upper:
        # 전하 과다 → 산성 잔기(D, E) 도입으로 순전하 감소
        diagnosis = (
            "Net charge exceeds +0.8 safety threshold. Excessive positive charge drives "
            "megalin-mediated renal tubular reabsorption, concentrating radioactivity in kidneys."
        )
        strategy = (
            "Introduce acidic residues (D: Asp, E: Glu) at solvent-exposed variable positions "
            "to reduce net positive charge below Dotatate reference (+0.707). Targeting positions "
            "outside the receptor-binding interface to avoid disrupting SSTR2 contacts."
        )
        for pos in mutable:
            if seq_list[pos] not in "DE" and rng.random() < 0.6:
                old_aa = seq_list[pos]
                new_aa = "D" if rng.random() < 0.5 else "E"
                seq_list[pos] = new_aa
                mutations.append(
                    f"Position {pos+1}: {old_aa} → {new_aa} "
                    f"(acidic — reduces net charge by ~1.0 unit at pH 7.4)"
                )
                if len(mutations) >= 2:
                    break

    elif "COMPOSITE" in reason_upper or "CSS" in reason_upper or "STABILITY" in reason_upper:
        # CSS 낮음 → 프로테아제 저항성/응집 개선
        diagnosis = (
            "Composite Stability Score below 0.3. Multiple PK indicators are suboptimal: "
            "likely high protease susceptibility (chymotrypsin-like cleavage at F/W/Y/L) "
            "and/or elevated aggregation tendency."
        )
        strategy = (
            "Replace exposed hydrophobic residues (L, I, V, A) with polar alternatives (S, T, N) "
            "to reduce aggregation tendency (Boman index) while maintaining core binding affinity. "
            "Avoid creating new chymotrypsin (F/W/Y/L) or trypsin (K/R) cleavage sites."
        )
        polar_pool = list("STNQG")
        for pos in mutable:
            if seq_list[pos] in "LIVA" and rng.random() < 0.6:
                old_aa = seq_list[pos]
                new_aa = polar_pool[rng.randint(len(polar_pool))]
                seq_list[pos] = new_aa
                mutations.append(
                    f"Position {pos+1}: {old_aa} → {new_aa} "
                    f"(polar — reduces aggregation, improves solubility)"
                )
                if len(mutations) >= 2:
                    break
    else:
        # 기본 fallback: 일반적 최적화
        diagnosis = (
            "Rejected by review panel. Specific failure mode not clearly identified from "
            "the rejection notes — applying general optimization strategy."
        )
        strategy = (
            "Conservative single-point mutation at the most solvent-exposed variable position, "
            "substituting with a residue that improves the weakest predicted indicator."
        )
        if mutable:
            pos = mutable[rng.randint(len(mutable))]
            old_aa = seq_list[pos]
            new_aa = VARIABLE_AA_POOL[rng.randint(len(VARIABLE_AA_POOL))]
            while new_aa == old_aa:
                new_aa = VARIABLE_AA_POOL[rng.randint(len(VARIABLE_AA_POOL))]
            seq_list[pos] = new_aa
            mutations.append(
                f"Position {pos+1}: {old_aa} → {new_aa} (conservative optimization)"
            )

    # 변이가 없었으면 최소 1개 강제
    if not mutations and mutable:
        pos = mutable[rng.randint(len(mutable))]
        old_aa = seq_list[pos]
        new_aa = VARIABLE_AA_POOL[rng.randint(len(VARIABLE_AA_POOL))]
        while new_aa == old_aa:
            new_aa = VARIABLE_AA_POOL[rng.randint(len(VARIABLE_AA_POOL))]
        seq_list[pos] = new_aa
        mutations.append(f"Position {pos+1}: {old_aa} → {new_aa} (general optimization)")

    new_seq = "".join(seq_list)
    mutation_details = "\n".join(f"  - {m}" for m in mutations)

    return (
        f"[Design Generator — Redesign Analysis]\n\n"
        f"[1. REJECTION DIAGNOSIS]\n{diagnosis}\n\n"
        f"[2. DESIGN STRATEGY]\n{strategy}\n\n"
        f"[3. MUTATION DETAILS]\nOriginal: {seq}\n{mutation_details}\n\n"
        f"[4. NEW SEQUENCE]\n"
        f"Proposed: {new_seq} (preserving Cys ring + FWKT pharmacophore)\n"
        f"NEW_SEQUENCE: {new_seq}"
    )


def _mock_experiment_planner(seq, data):
    energy = data.get("SSTR2_BindingEnergy", "N/A")
    sel = data.get("Selectivity", "N/A")
    charge = data.get("NetCharge", "N/A")
    return (
        f"[Experiment Planner - Preclinical Protocol]\n"
        f"Candidate: {seq} | Binding: {energy} kcal/mol | Sel: {sel}x | Charge: {charge}\n\n"
        f"1. DOTA-NHS Coupling:\n"
        f"   - Conjugate DOTA-NHS ester to N-terminal alpha-amine\n"
        f"   - Conditions: pH 8.5 carbonate buffer, 25C, 2h\n"
        f"   - Purify by RP-HPLC (C18), confirm purity >95% by LC-MS\n\n"
        f"2. Lu-177 Radiolabeling:\n"
        f"   - DOTA-peptide + [Lu-177]LuCl3, 0.4M acetate buffer pH 5.0, 95C, 30min\n"
        f"   - Radiochemical purity >98% by iTLC-SG\n\n"
        f"3. In Vitro Binding:\n"
        f"   - AR42J cells (endogenous SSTR2): saturation binding for Kd\n"
        f"   - Selectivity panel: competition binding on CHO-SSTR1/3/4/5\n\n"
        f"4. In Vivo Biodistribution (AR42J xenograft nude mice, n=5/timepoint):\n"
        f"   - Timepoints: 1h, 4h, 24h, 48h, 72h post-injection\n"
        f"   - Measure %ID/g: tumor, kidneys, liver, blood, spleen\n"
        f"   - Calculate tumor-to-kidney ratio at each timepoint\n\n"
        f"5. Success Criteria:\n"
        f"   - Tumor uptake: >10 %ID/g at 4h\n"
        f"   - 72h retention: tumor %ID/g decrease <50% vs 4h\n"
        f"   - Tumor/kidney ratio: >1.2x vs [Lu-177]DOTATATE (20% improvement)"
    )


# =========================================================================
#  LangGraph 상태 정의 (State Schema)
# =========================================================================

class BoardroomState(TypedDict):
    # 전체 후보 리스트 (Top-100 from Module 02)
    candidates: list
    # 현재 평가 중인 후보 인덱스
    current_idx: int
    # 현재 평가 중인 서열 (설계 생성기가 변이시킬 수 있음)
    current_sequence: str
    # 현재 후보의 전체 데이터
    current_data: dict
    # MCP 공동 작업 노트 (에이전트 간 공유, 누적)
    mcp_notebook: Annotated[list, add]
    # 전체 토론 기록 (Transcript, 누적)
    transcript: Annotated[list, add]
    # 구조 분석가 판정 (PASS / REJECT)
    structure_verdict: str
    # PK 전문가 판정 (PASS / REJECT)
    pk_verdict: str
    # 현재 후보의 재설계 반복 횟수
    iteration: int
    # 최대 재설계 횟수
    max_iterations: int
    # 최종 승인 후보 (누적)
    approved: Annotated[list, add]
    # 최종 보고서 파트 (누적)
    report_parts: Annotated[list, add]


# =========================================================================
#  에이전트 노드 함수 (Agent Node Functions)
# =========================================================================

# 글로벌 LLM 인스턴스 (main에서 초기화)
_LLM = None


def _format_candidate_context(state):
    """현재 후보의 컨텍스트 문자열을 생성한다."""
    seq = state["current_sequence"]
    data = state["current_data"]
    # MCP 노트에서 최근 3개만 표시 (컨텍스트 길이 관리)
    recent_notes = state.get("mcp_notebook", [])[-3:]
    notes_text = "\n".join(
        f"  [{n.get('agent', '?')}] {n.get('note', '')[:200]}"
        for n in recent_notes
    ) if recent_notes else "  (아직 없음)"

    return (
        f"=== Current Candidate ===\n"
        f"Sequence: {seq}\n"
        f"Length: {len(seq)} AA\n"
        f"SSTR2 Binding Energy: {data.get('SSTR2_BindingEnergy', 'N/A')} kcal/mol\n"
        f"Selectivity: {data.get('Selectivity', 'N/A')}x (vs non-target)\n"
        f"Net Charge (pH 7.4): {data.get('NetCharge', 'N/A')}\n"
        f"pI: {data.get('pI', 'N/A')}\n"
        f"Hydrophobicity: {data.get('Hydrophobicity', 'N/A')}\n"
        f"\n=== Module 01: 5-Indicator Stability Profile ===\n"
        f"Protease Cleavage Sites: {data.get('ProteaseSites', 'N/A')}\n"
        f"Aggregation Score: {data.get('AggregationScore', 'N/A')}\n"
        f"PPB Score (LogP/TPSA): {data.get('PPBScore', 'N/A')}\n"
        f"Rigidity (Rotatable Bonds): {data.get('Rigidity', 'N/A')}\n"
        f"Composite Stability Score: {data.get('CompositeStabilityScore', 'N/A')}\n"
        f"\n=== MCP Shared Notebook (recent) ===\n{notes_text}"
    )


def evidence_curator_node(state: dict) -> dict:
    """
    [증거 큐레이터] 논문/PDB 기반 지식을 제공한다.
    Trp-Lys 마스터키 메커니즘이 서열에서 잘 구현되었는지 분석.
    """
    seq = state["current_sequence"]
    data = state["current_data"]
    idx = state.get("current_idx", 0)
    print(f"    [증거 큐레이터] 후보 #{idx+1} {seq} 분석 중...", flush=True)

    context = _format_candidate_context(state)
    prompt = f"{context}\n\n위 서열을 분석하여 SSTR2 결합에 유리한 구조적 특성과 관련 문헌 레퍼런스를 제공하라."

    response = _call_llm(_LLM, PROMPT_EVIDENCE_CURATOR, prompt)
    if response is None:
        response = _mock_evidence_curator(seq, data)

    print(f"    [증거 큐레이터] 완료", flush=True)
    return {
        "mcp_notebook": [{"agent": "Evidence Curator", "note": response[:500]}],
        "transcript": [{"agent": "Evidence Curator", "sequence": seq, "message": response}],
    }


def structure_analyst_node(state: dict) -> dict:
    """
    [구조 분석가] 성공 지표 1(100배 정밀타격) + 2(72시간 생존력)를 검증한다.
    VERDICT: PASS 또는 REJECT를 명시적으로 반환.
    """
    seq = state["current_sequence"]
    data = state["current_data"]
    print(f"    [구조 분석가] {seq} 평가 중... (선택성/결합력)", flush=True)

    context = _format_candidate_context(state)
    prompt = (
        f"{context}\n\n"
        f"위 서열을 [성공 지표 1: 100배 정밀타격]과 [성공 지표 2: 72시간 생존력] "
        f"관점에서 엄격히 평가하라. 반드시 VERDICT: PASS 또는 VERDICT: REJECT로 결론을 내려라."
    )

    response = _call_llm(_LLM, PROMPT_STRUCTURE_ANALYST, prompt)
    if response is None:
        response = _mock_structure_analyst(seq, data)

    # 판정 파싱
    verdict = "REJECT"
    if "VERDICT: PASS" in response.upper() or "VERDICT:PASS" in response.upper():
        verdict = "PASS"

    print(f"    [구조 분석가] → {verdict}", flush=True)
    return {
        "structure_verdict": verdict,
        "mcp_notebook": [{"agent": "Structure Analyst", "note": f"[{verdict}] {response[:400]}"}],
        "transcript": [{"agent": "Structure Analyst", "sequence": seq, "verdict": verdict, "message": response}],
    }


def pk_expert_node(state: dict) -> dict:
    """
    [PK 전문가] 성공 지표 3(신장 독성 20% 개선)을 검증한다.
    VERDICT: PASS 또는 REJECT를 명시적으로 반환.
    """
    seq = state["current_sequence"]
    data = state["current_data"]
    print(f"    [PK 전문가] {seq} 평가 중... (전하/안정성/신장독성)", flush=True)

    context = _format_candidate_context(state)
    prompt = (
        f"{context}\n\n"
        f"Evaluate this candidate's PK profile using ALL 5 stability indicators "
        f"(Protease, Aggregation, Kidney/NetCharge, PPB, Rigidity) and the "
        f"Composite Stability Score. Compare NetCharge against Dotatate reference (+0.707). "
        f"Conclude with VERDICT: PASS or VERDICT: REJECT."
    )

    response = _call_llm(_LLM, PROMPT_PK_EXPERT, prompt)
    if response is None:
        response = _mock_pk_expert(seq, data)

    verdict = "REJECT"
    if "VERDICT: PASS" in response.upper() or "VERDICT:PASS" in response.upper():
        verdict = "PASS"

    print(f"    [PK 전문가] → {verdict}", flush=True)
    return {
        "pk_verdict": verdict,
        "mcp_notebook": [{"agent": "PK Expert", "note": f"[{verdict}] {response[:400]}"}],
        "transcript": [{"agent": "PK Expert", "sequence": seq, "verdict": verdict, "message": response}],
    }


def design_generator_node(state: dict) -> dict:
    """
    [설계 생성기] 반려된 서열의 실패 사유를 분석하고,
    과학적 근거에 기반한 변이 전략을 제시한 뒤 새로운 서열을 제안한다.
    핵심 뼈대(Cys 고리 + FWKT)는 절대 변경하지 않는다.
    """
    seq = state["current_sequence"]
    data = state["current_data"]
    iteration = state["iteration"]
    print(f"    [설계 생성기] {seq} 재설계 중... (변이 #{iteration+1})", flush=True)

    # 반려 사유 수집 — 모든 이전 에이전트 노트를 포함
    all_notes = state.get("mcp_notebook", [])
    rejection_reasons = "\n".join(
        f"  [{n.get('agent', '?')}] {n.get('note', '')}"
        for n in all_notes if "REJECT" in n.get("note", "")
    )
    # 증거 큐레이터의 구조 분석도 포함 (참고용)
    curator_notes = "\n".join(
        f"  [{n.get('agent', '?')}] {n.get('note', '')[:300]}"
        for n in all_notes if n.get("agent") == "Evidence Curator"
    )

    context = _format_candidate_context(state)
    prompt = (
        f"{context}\n\n"
        f"=== Rejection History (Iteration {iteration + 1}) ===\n"
        f"{rejection_reasons}\n\n"
        f"=== Structural Evidence (from Evidence Curator) ===\n"
        f"{curator_notes}\n\n"
        f"Based on the rejection reasons above, diagnose the specific failure mode, "
        f"explain your redesign strategy with SAR rationale, list each mutation with justification, "
        f"and propose a new sequence. Follow the structured format in your system prompt."
    )

    response = _call_llm(_LLM, PROMPT_DESIGN_GENERATOR, prompt)
    if response is None:
        response = _mock_design_generator(seq, rejection_reasons)

    # 새 서열 파싱
    new_seq = seq  # 파싱 실패 시 원본 유지
    for line in response.split("\n"):
        if "NEW_SEQUENCE:" in line.upper():
            candidate = line.split(":", 1)[1].strip().upper()
            # 유효성 검증: FWKT 포함, C 2개 이상
            if "FWKT" in candidate and candidate.count("C") >= 2:
                # K, R이 FWKT 외부에 없는지 확인
                temp = candidate.replace("FWKT", "XXXX")
                if "K" not in temp and "R" not in temp:
                    new_seq = candidate
            break

    # MCP 노트에 핵심 전략 요약 포함 (다른 에이전트가 참고)
    strategy_summary = ""
    for line in response.split("\n"):
        if "DESIGN STRATEGY" in line.upper() or "STRATEGY" in line.upper():
            # 전략 섹션의 다음 몇 줄을 가져옴
            idx = response.find(line)
            strategy_summary = response[idx:idx+300].strip()
            break
    if not strategy_summary:
        strategy_summary = f"변이: {seq} -> {new_seq}"

    print(f"    [설계 생성기] → {seq} → {new_seq}", flush=True)
    return {
        "current_sequence": new_seq,
        "iteration": iteration + 1,
        "mcp_notebook": [{"agent": "Design Generator",
                          "note": f"[Iter {iteration+1}] {seq} → {new_seq} | {strategy_summary[:200]}"}],
        "transcript": [{"agent": "Design Generator", "sequence": new_seq,
                        "iteration": iteration + 1, "message": response}],
    }


def experiment_planner_node(state: dict) -> dict:
    """
    [실험 계획가] 최종 승인된 서열에 대한 실험 계획서를 작성한다.
    DOTA 부착 전략, Lu-177 표지, In vivo 검증 프로토콜 포함.
    """
    seq = state["current_sequence"]
    data = state["current_data"]
    print(f"    [실험 계획가] {seq} APPROVED! 실험 계획서 작성 중...", flush=True)
    context = _format_candidate_context(state)

    prompt = (
        f"{context}\n\n"
        f"이 서열은 구조 분석가와 PK 전문가의 심사를 통과했습니다. "
        f"DOTA 부착, Lu-177 표지, In vivo 실험 계획을 작성하라."
    )

    response = _call_llm(_LLM, PROMPT_EXPERIMENT_PLANNER, prompt)
    if response is None:
        response = _mock_experiment_planner(seq, data)

    approved_entry = {
        "sequence": seq,
        "data": data,
        "experiment_plan": response,
    }

    return {
        "approved": [approved_entry],
        "report_parts": [{"sequence": seq, "plan": response}],
        "mcp_notebook": [{"agent": "Experiment Planner", "note": f"실험 계획 완료: {seq}"}],
        "transcript": [{"agent": "Experiment Planner", "sequence": seq, "message": response}],
    }


def next_candidate_node(state: dict) -> dict:
    """다음 후보로 이동한다. 인덱스 증가, 상태 초기화."""
    new_idx = state["current_idx"] + 1
    total = len(state["candidates"])
    if new_idx < total:
        print(f"\n  --- 후보 #{new_idx+1}/{total}: "
              f"{state['candidates'][new_idx]['Sequence']} ---", flush=True)
    else:
        print(f"\n  === 모든 후보 평가 완료 ({total}개) ===", flush=True)
    candidates = state["candidates"]

    if new_idx < len(candidates):
        next_cand = candidates[new_idx]
        return {
            "current_idx": new_idx,
            "current_sequence": next_cand["Sequence"],
            "current_data": next_cand,
            "structure_verdict": "",
            "pk_verdict": "",
            "iteration": 0,
            "transcript": [{"agent": "System", "message": f"--- 후보 #{new_idx + 1} 평가 시작 ---"}],
        }
    else:
        return {
            "current_idx": new_idx,
            "transcript": [{"agent": "System", "message": "=== 모든 후보 평가 완료 ==="}],
        }


# =========================================================================
#  라우팅 함수 (Routing Functions)
# =========================================================================

def route_after_structure(state: dict) -> str:
    """구조 분석가 판정 후 라우팅."""
    if state["structure_verdict"] == "PASS":
        return "pk_expert"
    elif state["iteration"] >= state["max_iterations"]:
        # 최대 재설계 횟수 초과 → 이 후보 포기
        return "next_candidate"
    else:
        return "design_generator"


def route_after_pk(state: dict) -> str:
    """PK 전문가 판정 후 라우팅."""
    if state["pk_verdict"] == "PASS":
        return "experiment_planner"
    elif state["iteration"] >= state["max_iterations"]:
        return "next_candidate"
    else:
        return "design_generator"


def route_next_or_end(state: dict) -> str:
    """다음 후보가 있으면 계속, 없으면 종료."""
    if state["current_idx"] < len(state["candidates"]):
        return "evidence_curator"
    else:
        return END


# =========================================================================
#  그래프 구축 (Graph Construction)
# =========================================================================

def build_boardroom_graph():
    """
    LangGraph StateGraph로 5인 에이전트 순환 협업 그래프를 구축한다.

    [그래프 구조 - State Machine]
    START -> evidence_curator -> structure_analyst -> [분기]
      PASS  -> pk_expert -> [분기]
                PASS  -> experiment_planner -> next_candidate -> [분기]
                                                  더 있음 -> evidence_curator
                                                  끝     -> END
                REJECT -> design_generator -> structure_analyst (순환)
      REJECT -> design_generator -> structure_analyst (순환)
    max_iterations 초과 시 -> next_candidate (포기)
    """
    builder = StateGraph(BoardroomState)

    # 노드 등록
    builder.add_node("evidence_curator", evidence_curator_node)
    builder.add_node("structure_analyst", structure_analyst_node)
    builder.add_node("pk_expert", pk_expert_node)
    builder.add_node("design_generator", design_generator_node)
    builder.add_node("experiment_planner", experiment_planner_node)
    builder.add_node("next_candidate", next_candidate_node)

    # 엣지 연결
    builder.add_edge(START, "evidence_curator")
    builder.add_edge("evidence_curator", "structure_analyst")

    # 구조 분석가 -> 조건부 분기
    builder.add_conditional_edges(
        "structure_analyst", route_after_structure,
        {"pk_expert": "pk_expert", "design_generator": "design_generator",
         "next_candidate": "next_candidate"},
    )
    # PK 전문가 -> 조건부 분기
    builder.add_conditional_edges(
        "pk_expert", route_after_pk,
        {"experiment_planner": "experiment_planner", "design_generator": "design_generator",
         "next_candidate": "next_candidate"},
    )
    # 설계 생성기 -> 구조 분석가 (재평가 순환)
    builder.add_edge("design_generator", "structure_analyst")
    # 실험 계획가 -> 다음 후보
    builder.add_edge("experiment_planner", "next_candidate")
    # 다음 후보 -> 조건부 분기 (계속 or 종료)
    builder.add_conditional_edges(
        "next_candidate", route_next_or_end,
        {"evidence_curator": "evidence_curator", END: END},
    )

    return builder.compile()


# =========================================================================
#  Fallback: LangGraph 없이 순차 실행
# =========================================================================

def _call_design_generator(current_seq, data, rejecting_agent, rejection_resp,
                           evidence_resp, iteration, transcript):
    """
    Design Generator를 호출하여 반려 사유 분석 + 변이 서열을 생성한다.
    LLM과 Mock 모두 동일한 풍부한 컨텍스트를 전달한다.
    """
    # Mock 응답 생성 (항상 — LLM 실패 시 fallback)
    dg_resp = _mock_design_generator(current_seq, rejection_resp)

    # LLM 호출 (가능하면)
    if _LLM:
        llm_ctx = (
            f"=== Candidate Data ===\n"
            f"Sequence: {current_seq}\n"
            f"Binding Energy: {data.get('SSTR2_BindingEnergy', 'N/A')} kcal/mol\n"
            f"Selectivity: {data.get('Selectivity', 'N/A')}x\n"
            f"Net Charge: {data.get('NetCharge', 'N/A')}\n"
            f"Composite Stability Score: {data.get('CompositeStabilityScore', 'N/A')}\n"
            f"Protease Sites: {data.get('ProteaseSites', 'N/A')}\n"
            f"Aggregation: {data.get('AggregationScore', 'N/A')}\n"
            f"PPB: {data.get('PPBScore', 'N/A')}\n"
            f"Rigidity: {data.get('Rigidity', 'N/A')}\n\n"
            f"=== Evidence Curator's Analysis ===\n"
            f"{evidence_resp[:500]}\n\n"
            f"=== Rejection by {rejecting_agent} (Iteration {iteration + 1}) ===\n"
            f"{rejection_resp}\n\n"
            f"Diagnose the failure, explain your redesign strategy with SAR rationale, "
            f"and propose a new sequence following the structured format."
        )
        real_resp = _call_llm(_LLM, PROMPT_DESIGN_GENERATOR, llm_ctx)
        if real_resp:
            dg_resp = real_resp

    # 서열 파싱
    new_seq = current_seq
    for line in dg_resp.split("\n"):
        if "NEW_SEQUENCE:" in line.upper():
            candidate_seq = line.split(":", 1)[1].strip().upper()
            if "FWKT" in candidate_seq and candidate_seq.count("C") >= 2:
                temp = candidate_seq.replace("FWKT", "XXXX")
                if "K" not in temp and "R" not in temp:
                    new_seq = candidate_seq
            break

    # 콘솔 출력 — 핵심 전략만 요약
    print(f"    설계 생성기 (변이 #{iteration + 1}):")
    print(f"      {current_seq} → {new_seq}")
    for line in dg_resp.split("\n"):
        stripped = line.strip()
        if stripped.startswith("Position") or stripped.startswith("- Position"):
            print(f"      {stripped[:80]}")

    transcript.append({"agent": "Design Generator", "sequence": new_seq,
                       "iteration": iteration + 1, "message": dg_resp})
    return new_seq, dg_resp


def run_sequential_boardroom(candidates):
    """
    LangGraph가 없을 때 순차적으로 에이전트를 실행한다.
    동일한 로직을 for 루프로 구현.
    """
    transcript = []
    approved = []
    report_parts = []

    for c_idx, cand in enumerate(candidates):
        seq = cand["Sequence"]
        data = cand
        print(f"\n  --- 후보 #{c_idx + 1}/{len(candidates)}: {seq} ---")
        transcript.append({"agent": "System", "message": f"--- 후보 #{c_idx + 1} ---"})

        iteration = 0
        current_seq = seq
        evidence_resp = ""  # 증거 큐레이터 응답 보존 (설계 생성기 참고용)

        while iteration <= MAX_REDESIGN_ITERATIONS:
            # 1) 증거 큐레이터
            ec_resp = _mock_evidence_curator(current_seq, data)
            transcript.append({"agent": "Evidence Curator", "sequence": current_seq, "message": ec_resp})
            if _LLM:
                ctx = (
                    f"Sequence: {current_seq}\n"
                    f"Binding Energy: {data.get('SSTR2_BindingEnergy')} kcal/mol\n"
                    f"Selectivity: {data.get('Selectivity')}x\n"
                    f"Net Charge: {data.get('NetCharge')}\n"
                    f"Analyze structural compatibility with SSTR2 binding pocket (PDB 7T10)."
                )
                real_resp = _call_llm(_LLM, PROMPT_EVIDENCE_CURATOR, ctx)
                if real_resp:
                    ec_resp = real_resp
                    transcript[-1]["message"] = real_resp
            evidence_resp = ec_resp

            # 2) 구조 분석가
            sa_resp = _mock_structure_analyst(current_seq, data)
            if _LLM:
                ctx = (
                    f"Sequence: {current_seq}\n"
                    f"Binding Energy: {data.get('SSTR2_BindingEnergy')} kcal/mol\n"
                    f"Selectivity: {data.get('Selectivity')}x\n"
                    f"Rigidity: {data.get('Rigidity')} rotatable bonds\n"
                    f"Evidence Curator notes: {ec_resp[:300]}\n\n"
                    f"Evaluate Criterion 1 (100x Selectivity) and Criterion 2 (72h Retention). "
                    f"Conclude with VERDICT: PASS or VERDICT: REJECT."
                )
                real_resp = _call_llm(_LLM, PROMPT_STRUCTURE_ANALYST, ctx)
                if real_resp:
                    sa_resp = real_resp
            s_verdict = "PASS" if "VERDICT: PASS" in sa_resp.upper() else "REJECT"
            print(f"    구조 분석가: {s_verdict}")
            transcript.append({"agent": "Structure Analyst", "sequence": current_seq,
                             "verdict": s_verdict, "message": sa_resp})

            if s_verdict == "REJECT":
                if iteration >= MAX_REDESIGN_ITERATIONS:
                    print(f"    -> 최대 재설계 초과, 포기")
                    break
                current_seq, _ = _call_design_generator(
                    current_seq, data, "Structure Analyst", sa_resp,
                    evidence_resp, iteration, transcript)
                iteration += 1
                continue

            # 3) PK 전문가
            pk_resp = _mock_pk_expert(current_seq, data)
            if _LLM:
                ctx = (
                    f"Sequence: {current_seq}\n"
                    f"Net Charge: {data.get('NetCharge')}\n"
                    f"pI: {data.get('pI')}\n"
                    f"Composite Stability Score: {data.get('CompositeStabilityScore')}\n"
                    f"Protease Sites: {data.get('ProteaseSites')}\n"
                    f"Aggregation Score: {data.get('AggregationScore')}\n"
                    f"PPB Score: {data.get('PPBScore')}\n"
                    f"Rigidity: {data.get('Rigidity')}\n\n"
                    f"Evaluate all 5 PK indicators and Criterion 3 (Kidney Safety). "
                    f"Conclude with VERDICT: PASS or VERDICT: REJECT."
                )
                real_resp = _call_llm(_LLM, PROMPT_PK_EXPERT, ctx)
                if real_resp:
                    pk_resp = real_resp
            p_verdict = "PASS" if "VERDICT: PASS" in pk_resp.upper() else "REJECT"
            print(f"    PK 전문가: {p_verdict}")
            transcript.append({"agent": "PK Expert", "sequence": current_seq,
                             "verdict": p_verdict, "message": pk_resp})

            if p_verdict == "REJECT":
                if iteration >= MAX_REDESIGN_ITERATIONS:
                    print(f"    -> 최대 재설계 초과, 포기")
                    break
                current_seq, _ = _call_design_generator(
                    current_seq, data, "PK Expert", pk_resp,
                    evidence_resp, iteration, transcript)
                iteration += 1
                continue

            # 4) 실험 계획가 (둘 다 PASS)
            ep_resp = _mock_experiment_planner(current_seq, data)
            if _LLM:
                real_resp = _call_llm(_LLM, PROMPT_EXPERIMENT_PLANNER,
                    f"Approved sequence: {current_seq}\n"
                    f"Binding: {data.get('SSTR2_BindingEnergy')} kcal/mol\n"
                    f"Selectivity: {data.get('Selectivity')}x\n"
                    f"Charge: {data.get('NetCharge')}\n"
                    f"Write a preclinical experiment protocol.")
                if real_resp:
                    ep_resp = real_resp
            transcript.append({"agent": "Experiment Planner", "sequence": current_seq,
                             "message": ep_resp})
            approved.append({"sequence": current_seq, "data": data, "experiment_plan": ep_resp})
            report_parts.append({"sequence": current_seq, "plan": ep_resp})
            print(f"    ==> APPROVED!")
            break

    return transcript, approved, report_parts


# =========================================================================
#  보고서 생성 (Report Generation)
# =========================================================================

def save_transcript(transcript, path):
    """Save agent discussion transcript as a text log."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("  SSTR2 AI Scientist - Boardroom Transcript\n")
        f.write(f"  Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 70 + "\n\n")
        for entry in transcript:
            agent = entry.get("agent", "Unknown")
            seq = entry.get("sequence", "")
            verdict = entry.get("verdict", "")
            msg = entry.get("message", "")
            f.write(f"[{agent}]")
            if seq:
                f.write(f" (Seq: {seq})")
            if verdict:
                f.write(f" [{verdict}]")
            f.write("\n")
            for line in msg.split("\n"):
                f.write(f"  {line}\n")
            f.write("\n" + "-" * 50 + "\n\n")


def _generate_selection_rationale(approved_entry):
    """Generate a structured selection rationale for an approved candidate."""
    seq = approved_entry["sequence"]
    d = approved_entry["data"]
    energy = d.get("SSTR2_BindingEnergy", "N/A")
    sel = d.get("Selectivity", "N/A")
    charge = d.get("NetCharge", "N/A")
    pi_val = d.get("pI", "N/A")
    hydro = d.get("Hydrophobicity", "N/A")
    length = d.get("Length", len(seq))
    composite = d.get("CompositeScore", "N/A")

    # Structural analysis
    c_pos = [i for i, a in enumerate(seq) if a == "C"]
    ring_info = f"Cys{c_pos[0]+1}-Cys{c_pos[1]+1}" if len(c_pos) >= 2 else "N/A"
    ring_size = (c_pos[1] - c_pos[0] - 1) if len(c_pos) >= 2 else 0
    fwkt_pos = seq.find("FWKT")

    # Charge assessment
    if isinstance(charge, (int, float)):
        if charge < 0:
            charge_assessment = f"Net charge ({charge:+.4f}) is negative — substantially improved over Dotatate (+0.707). Reduced megalin-mediated renal reabsorption expected."
        elif charge <= 0.707:
            charge_assessment = f"Net charge ({charge:+.4f}) is lower than Dotatate (+0.707) — kidney safety improved."
        else:
            charge_assessment = f"Net charge ({charge:+.4f}) is comparable to Dotatate (+0.707) — standard lysine co-infusion protocol applies."
    else:
        charge_assessment = "Charge data not available."

    return (
        f"**Structural Basis:**\n"
        f"- Scaffold: {length} AA, disulfide ring {ring_info} (ring size: {ring_size} residues)\n"
        f"- FWKT pharmacophore at position {fwkt_pos+1}-{fwkt_pos+4}: Trp inserts into SSTR2 TM3-TM7 "
        f"hydrophobic pocket (PDB 7T10); Lys forms salt bridge with Asp3.32\n"
        f"- Beta-turn constraint from disulfide ring reduces conformational entropy penalty upon binding\n\n"
        f"**Criterion 1 — 100x Selectivity:** Selectivity = {sel}x (threshold: >= 100x). "
        f"The FWKT motif provides SSTR2-specific interactions that are poorly accommodated by "
        f"SSTR1 (narrow Trp pocket), SSTR3 (ECL2 steric clash), SSTR4 (divergent binding site), "
        f"and SSTR5 (weakened Lys salt bridge).\n\n"
        f"**Criterion 2 — 72h Retention:** Binding energy = {energy} kcal/mol (threshold: <= -7.0). "
        f"Strong binding yields k_off < 3.86e-6 /s, enabling >72h target residence for effective "
        f"Lu-177 radiation delivery.\n\n"
        f"**Criterion 3 — Kidney Safety:** {charge_assessment}\n\n"
        f"**Composite Score:** {composite} (binding energy - log10(selectivity); lower = better)"
    )


def save_report(approved, report_parts, transcript, path):
    """Generate the final winning peptides report in English Markdown."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("# SSTR2 AI Scientist — Final Winning Peptides Report\n\n")
        f.write(f"**Generated:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("## 1. Project Overview\n\n")
        f.write("This report presents the final candidates selected by the 5-agent AI Boardroom system "
                "for next-generation SSTR2-targeted Peptide Receptor Radionuclide Therapy (PRRT). "
                "Each candidate has been evaluated against three mandatory success criteria by "
                "specialized AI agents (Evidence Curator, Structure Analyst, PK Expert, "
                "Design Generator, Experiment Planner).\n\n")

        f.write("## 2. Success Criteria\n\n")
        f.write("| # | Criterion | Threshold | Rationale |\n")
        f.write("|---|-----------|-----------|----------|\n")
        f.write("| 1 | **Target Selectivity** | Selectivity >= 100x vs SSTR1/3/4/5 | "
                "Minimize off-target radiation to healthy tissue |\n")
        f.write("| 2 | **72h Retention** | Binding energy <= -7.0 kcal/mol | "
                "k_off < 3.86e-6 /s for sustained Lu-177 irradiation |\n")
        f.write("| 3 | **Kidney Safety** | Net charge <= Dotatate (+0.707) | "
                "Reduce megalin-mediated renal tubular reabsorption |\n\n")

        f.write(f"## 3. Approved Candidates ({len(approved)})\n\n")
        if approved:
            f.write("| # | Sequence | Binding Energy (kcal/mol) | Selectivity | "
                    "Net Charge | Composite Score |\n")
            f.write("|---|----------|--------------------------|-------------|------------|----------------|\n")
            for i, a in enumerate(approved, 1):
                d = a["data"]
                f.write(
                    f"| {i} | `{a['sequence']}` | "
                    f"{d.get('SSTR2_BindingEnergy', 'N/A')} | "
                    f"{d.get('Selectivity', 'N/A')}x | "
                    f"{d.get('NetCharge', 'N/A')} | "
                    f"{d.get('CompositeScore', 'N/A')} |\n"
                )
            f.write("\n")

            # --- Selection Rationale (NEW) ---
            f.write("## 4. Selection Rationale\n\n")
            f.write("Each candidate was approved based on the following structural, biophysical, "
                    "and pharmacokinetic evidence:\n\n")
            for i, a in enumerate(approved, 1):
                f.write(f"### Candidate #{i}: `{a['sequence']}`\n\n")
                rationale = _generate_selection_rationale(a)
                f.write(rationale)
                f.write("\n\n---\n\n")

            # --- Experiment Plans ---
            f.write("## 5. Preclinical Experiment Plans\n\n")
            for i, rp in enumerate(report_parts, 1):
                f.write(f"### Candidate #{i}: `{rp['sequence']}`\n\n")
                f.write("```\n")
                f.write(rp["plan"])
                f.write("\n```\n\n")
        else:
            f.write("*No candidates were approved. Consider adjusting filter criteria "
                    "or evaluating more candidates.*\n\n")
            f.write("## 4. Selection Rationale\n\nN/A\n\n")
            f.write("## 5. Preclinical Experiment Plans\n\nN/A\n\n")

        # --- Agent Statistics ---
        f.write("## 6. Agent Evaluation Summary\n\n")
        agent_stats = {}
        for entry in transcript:
            agent = entry.get("agent", "")
            verdict = entry.get("verdict", "")
            if verdict:
                if agent not in agent_stats:
                    agent_stats[agent] = {"PASS": 0, "REJECT": 0}
                agent_stats[agent][verdict] = agent_stats[agent].get(verdict, 0) + 1
        if agent_stats:
            f.write("| Agent | PASS | REJECT |\n")
            f.write("|-------|------|--------|\n")
            for agent, stats in agent_stats.items():
                f.write(f"| {agent} | {stats.get('PASS', 0)} | {stats.get('REJECT', 0)} |\n")
            f.write("\n")

        f.write("---\n")
        f.write("*Generated by SSTR2 AI Scientist Boardroom System*\n")


# =========================================================================
#  메인 실행부 (Main Execution)
# =========================================================================

def main():
    global _LLM

    print("=" * 70)
    print("  SSTR2 AI Scientist - Module 03: 5인 AI 드림팀 원탁회의")
    print("=" * 70)

    start_time = time.time()

    # ---- Step 1: 데이터 로드 ----
    input_path = OUTPUT_DIR / "SSTR2_Top100_Boardroom_Ready.csv"
    transcript_path = OUTPUT_DIR / "Boardroom_Transcript.log"
    report_path = OUTPUT_DIR / "Final_Winning_Peptides_Report.md"

    if not input_path.exists():
        print(f"  [오류] 입력 파일 없음: {input_path}")
        print(f"  -> 먼저 02_train_surrogate.py를 실행하세요.")
        sys.exit(1)

    print(f"\n>>> Step 1: Top 후보 로드")
    df = pd.read_csv(input_path, encoding="utf-8")
    # 상위 N개만 원탁회의에 투입 (LLM 호출 비용 관리)
    top_n = min(BOARDROOM_TOP_N, len(df))
    candidates = df.head(top_n).to_dict("records")
    print(f"  로드: {len(df)}개 중 상위 {top_n}개를 원탁회의에 투입")

    # ---- Step 2: LLM 초기화 ----
    print(f"\n>>> Step 2: LLM 백엔드 초기화")
    print(f"  Ollama 모델: {OLLAMA_MODEL}")
    print(f"  Ollama URL: {OLLAMA_BASE_URL}")
    _LLM = _create_llm()

    # ---- Step 3: 원탁회의 실행 ----
    print(f"\n>>> Step 3: 원탁회의 시작 (후보 {top_n}개, 최대 재설계 {MAX_REDESIGN_ITERATIONS}회)")
    print(f"  에이전트: 증거큐레이터, 구조분석가, PK전문가, 설계생성기, 실험계획가")

    if LANGGRAPH_OK:
        # LangGraph State Machine 모드
        print(f"  [모드] LangGraph StateGraph (순환 구조)")
        graph = build_boardroom_graph()

        initial_state = {
            "candidates": candidates,
            "current_idx": 0,
            "current_sequence": candidates[0]["Sequence"],
            "current_data": candidates[0],
            "mcp_notebook": [],
            "transcript": [{"agent": "System", "message": "=== 원탁회의 개시 ==="}],
            "structure_verdict": "",
            "pk_verdict": "",
            "iteration": 0,
            "max_iterations": MAX_REDESIGN_ITERATIONS,
            "approved": [],
            "report_parts": [],
        }

        print(f"\n  LangGraph 실행 중...")
        # recursion_limit: 후보 N개 x (최대 11노드/후보) + 여유분
        rec_limit = top_n * 12 + 10
        try:
            final_state = graph.invoke(
                initial_state,
                {"recursion_limit": rec_limit},
            )
            transcript = final_state.get("transcript", [])
            approved = final_state.get("approved", [])
            report_parts = final_state.get("report_parts", [])
        except Exception as e:
            print(f"  [경고] LangGraph 실행 오류: {e}")
            print(f"  -> 순차 모드로 전환합니다.")
            transcript, approved, report_parts = run_sequential_boardroom(candidates)
    else:
        # Fallback: 순차 실행 모드
        print(f"  [모드] 순차 실행 (langgraph 미설치)")
        transcript, approved, report_parts = run_sequential_boardroom(candidates)

    # ---- Step 4: 결과 저장 ----
    print(f"\n>>> Step 4: 결과 저장")

    save_transcript(transcript, transcript_path)
    print(f"  토론 기록: {transcript_path}")

    save_report(approved, report_parts, transcript, report_path)
    print(f"  최종 보고서: {report_path}")

    elapsed = time.time() - start_time

    # ---- 최종 요약 ----
    print(f"\n{'=' * 70}")
    print(f"  >>> 원탁회의 완료!")
    print(f"  평가 후보: {top_n}개")
    print(f"  최종 승인: {len(approved)}개")
    print(f"  소요 시간: {elapsed:.2f}초")
    print(f"{'=' * 70}")

    if approved:
        print(f"\n  [최종 승인 위닝 펩타이드]")
        for i, a in enumerate(approved, 1):
            d = a["data"]
            print(
                f"  #{i} {a['sequence']} | "
                f"결합E={d.get('SSTR2_BindingEnergy', '?')} | "
                f"선택성={d.get('Selectivity', '?')}x | "
                f"전하={d.get('NetCharge', '?')}"
            )
    else:
        print(f"\n  [경고] 승인된 후보가 없습니다.")
        print(f"  -> BOARDROOM_TOP_N을 늘리거나 필터 기준을 조정하세요.")

    print(f"\n  [출력 파일]")
    print(f"  - {transcript_path}")
    print(f"  - {report_path}")

    return approved


if __name__ == "__main__":
    main()
