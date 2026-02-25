#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  SSTR2 AI Scientist - Real-time Mission Control Dashboard (dashboard.py)
==========================================================================

  Streamlit 기반 실시간 대시보드.
  파이프라인 진행 상태, 에이전트 활동, 후보 랭킹을 브라우저에서 모니터링.

  [실행 방법]
    streamlit run dashboard.py

  [데이터 소스 — 기존 모듈 출력 파일을 읽기 전용으로 참조]
    output/SSTR2_Raw_Library.csv          (Module 00)
    output/SSTR2_Filtered_10K.csv         (Module 01)
    output/SSTR2_Top100_Boardroom_Ready.csv (Module 02)
    output/Boardroom_Transcript.log       (Module 03)
    output/Final_Winning_Peptides_Report.md (Module 03)
==========================================================================
"""

import re
import time
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path

# =========================================================================
#  상수 및 경로 (Constants & Paths)
# =========================================================================

OUTPUT_DIR = Path("output")

# 파이프라인 출력 파일 경로
FILE_RAW_LIBRARY = OUTPUT_DIR / "SSTR2_Raw_Library.csv"
FILE_FILTERED = OUTPUT_DIR / "SSTR2_Filtered_Top100.csv"
FILE_FILTERED_LEGACY = OUTPUT_DIR / "SSTR2_Filtered_10K.csv"
FILE_TOP100 = OUTPUT_DIR / "SSTR2_Top100_Boardroom_Ready.csv"
FILE_TRANSCRIPT = OUTPUT_DIR / "Boardroom_Transcript.log"
FILE_REPORT = OUTPUT_DIR / "Final_Winning_Peptides_Report.md"

# 3대 성공 지표 임계값 (02_train_surrogate.py 에서 정의한 것과 동일)
DOTATATE_DOCKING_SCORE = -8.5    # kcal/mol
MIN_RETENTION_ENERGY = -7.0      # kcal/mol (72시간 체류 최소 결합 에너지)
SELECTIVITY_THRESHOLD = 100.0    # 100배 이상 선택성
DOTATATE_CHARGE = 0.707          # Dotatate 순전하

# 에이전트 색상 매핑 (Agent color mapping)
AGENT_COLORS = {
    "Evidence Curator": "#3498db",   # 파란색
    "Structure Analyst": "#e67e22",  # 주황색
    "PK Expert": "#2ecc71",          # 초록색
    "Design Generator": "#9b59b6",   # 보라색
    "Experiment Planner": "#e74c3c", # 빨간색
    "System": "#95a5a6",             # 회색
}

# 에이전트 아이콘 매핑
AGENT_ICONS = {
    "Evidence Curator": "📚",
    "Structure Analyst": "🔬",
    "PK Expert": "💊",
    "Design Generator": "🧬",
    "Experiment Planner": "📋",
    "System": "⚙️",
}


# =========================================================================
#  데이터 로딩 함수 (Data Loading — 읽기 전용)
# =========================================================================

@st.cache_data(ttl=5)
def load_pipeline_status():
    """각 모듈의 출력 파일 존재 여부로 파이프라인 진행 상태를 판단한다."""
    steps = [
        {"name": "Module 00", "label": "Peptide Generation",
         "file": FILE_RAW_LIBRARY, "status": "waiting"},
        {"name": "Module 01", "label": "PK Filter",
         "file": FILE_FILTERED if FILE_FILTERED.exists() else FILE_FILTERED_LEGACY,
         "status": "waiting"},
        {"name": "Module 02", "label": "Surrogate Model",
         "file": FILE_TOP100, "status": "waiting"},
        {"name": "Module 03", "label": "AI Boardroom",
         "file": FILE_TRANSCRIPT, "status": "waiting"},
    ]
    # 마지막 완료된 단계의 다음 단계가 "진행중"
    last_done = -1
    for i, s in enumerate(steps):
        if s["file"].exists():
            s["status"] = "done"
            last_done = i

    # 아직 완료되지 않은 다음 단계를 "진행중"으로 표시
    if last_done < len(steps) - 1:
        steps[last_done + 1]["status"] = "running"

    return steps


@st.cache_data(ttl=5)
def load_raw_library():
    """Module 00 출력 전체 DataFrame을 로드한다."""
    if FILE_RAW_LIBRARY.exists():
        return pd.read_csv(FILE_RAW_LIBRARY, encoding="utf-8")
    return pd.DataFrame()


@st.cache_data(ttl=5)
def load_filtered():
    """Module 01 출력 전체 DataFrame을 로드한다."""
    if FILE_FILTERED.exists():
        return pd.read_csv(FILE_FILTERED, encoding="utf-8")
    if FILE_FILTERED_LEGACY.exists():
        return pd.read_csv(FILE_FILTERED_LEGACY, encoding="utf-8")
    return pd.DataFrame()


@st.cache_data(ttl=5)
def load_top100():
    """Module 02 출력(Top-100 후보)을 로드한다."""
    if FILE_TOP100.exists():
        return pd.read_csv(FILE_TOP100, encoding="utf-8")
    return pd.DataFrame()


@st.cache_data(ttl=5)
def load_transcript():
    """Boardroom_Transcript.log를 파싱하여 에이전트 활동 목록을 반환한다."""
    if not FILE_TRANSCRIPT.exists():
        return []

    text = FILE_TRANSCRIPT.read_text(encoding="utf-8")
    # 구분자 기준으로 블록 분리
    blocks = text.split("-" * 50)
    entries = []

    for block in blocks:
        block = block.strip()
        if not block:
            continue
        # [Agent] (Seq: XXX) [VERDICT] pattern parsing (supports both English and Korean)
        header_match = re.match(
            r"\[([^\]]+)\](?:\s*\((?:Seq|서열):\s*([^)]+)\))?(?:\s*\[([^\]]+)\])?",
            block
        )
        if header_match:
            agent = header_match.group(1)
            sequence = header_match.group(2) or ""
            verdict = header_match.group(3) or ""
            # 헤더 이후 메시지 본문
            msg_start = header_match.end()
            message = block[msg_start:].strip()
            entries.append({
                "agent": agent,
                "sequence": sequence,
                "verdict": verdict,
                "message": message,
            })

    return entries


@st.cache_data(ttl=5)
def load_report():
    """최종 보고서 마크다운을 로드한다."""
    if FILE_REPORT.exists():
        return FILE_REPORT.read_text(encoding="utf-8")
    return ""


# =========================================================================
#  유틸리티 함수 (Utility Functions)
# =========================================================================

def judge_metrics(row):
    """3대 성공 지표 통과 여부를 판정한다."""
    energy = row.get("SSTR2_BindingEnergy", 0)
    sel = row.get("Selectivity", 0)
    charge = row.get("NetCharge", 999)

    # 지표 1: 선택성 >= 100배
    metric1 = sel >= SELECTIVITY_THRESHOLD
    # 지표 2: 결합 에너지 <= -7.0 (72시간 체류)
    metric2 = energy <= MIN_RETENTION_ENERGY
    # 지표 3: 순전하 <= Dotatate (+0.707)
    metric3 = charge <= DOTATATE_CHARGE

    return metric1, metric2, metric3


def render_verdict_badge(passed):
    """PASS/FAIL 뱃지 HTML을 생성한다."""
    if passed:
        return '<span style="background:#2ecc71;color:white;padding:2px 8px;border-radius:4px;font-weight:bold;font-size:12px;">PASS</span>'
    else:
        return '<span style="background:#e74c3c;color:white;padding:2px 8px;border-radius:4px;font-weight:bold;font-size:12px;">FAIL</span>'


# =========================================================================
#  UI 컴포넌트: 상단 프로세스 바 (Pipeline Progress)
# =========================================================================

def render_pipeline_progress():
    """파이프라인 진행 상태를 상단 프로세스 바로 렌더링한다."""
    steps = load_pipeline_status()
    total = len(steps)
    done_count = sum(1 for s in steps if s["status"] == "done")

    # 전체 진행률
    progress_pct = done_count / total
    st.progress(progress_pct, text=f"Pipeline Progress: {done_count}/{total} steps completed")

    # 각 Step 상세 표시
    cols = st.columns(total)
    for i, (col, step) in enumerate(zip(cols, steps)):
        status = step["status"]
        with col:
            if status == "done":
                icon = "✅"
                color = "#2ecc71"
            elif status == "running":
                icon = "🔄"
                color = "#f39c12"
            else:
                icon = "⏳"
                color = "#95a5a6"

            st.markdown(
                f"""<div style="text-align:center;padding:8px;border-radius:8px;
                border:2px solid {color};background:rgba(0,0,0,0.02);">
                <div style="font-size:24px;">{icon}</div>
                <div style="font-weight:bold;color:{color};font-size:13px;">{step['name']}</div>
                <div style="font-size:11px;color:#666;">{step['label']}</div>
                </div>""",
                unsafe_allow_html=True,
            )

    # 추가 통계
    df_raw = load_raw_library()
    df_filt = load_filtered()
    df_top = load_top100()

    m1, m2, m3 = st.columns(3)
    with m1:
        st.metric("Generated Sequences", f"{len(df_raw):,}" if not df_raw.empty else "—")
    with m2:
        st.metric("PK-Filtered", f"{len(df_filt):,}" if not df_filt.empty else "—")
    with m3:
        st.metric("Top Candidates", f"{len(df_top):,}" if not df_top.empty else "—")


# =========================================================================
#  UI 컴포넌트: 사이드바 Agent Monitor
# =========================================================================

def render_agent_monitor():
    """사이드바에 에이전트 활동 로그를 실시간 스트리밍한다."""
    st.sidebar.markdown("## 🤖 Agent Monitor")

    entries = load_transcript()

    if not entries:
        st.sidebar.info("Boardroom 로그 대기 중...\n\n03_ai_boardroom.py 실행 후 로그가 표시됩니다.")
        return

    # 에이전트별 최신 상태 요약
    st.sidebar.markdown("### Agent Status")
    # 각 에이전트의 마지막 활동 수집
    latest = {}
    for e in entries:
        agent = e["agent"]
        if agent != "System":
            latest[agent] = e

    for agent, entry in latest.items():
        icon = AGENT_ICONS.get(agent, "🔹")
        color = AGENT_COLORS.get(agent, "#666")
        verdict = entry.get("verdict", "")
        seq = entry.get("sequence", "")

        verdict_html = ""
        if verdict == "PASS":
            verdict_html = ' <span style="color:#2ecc71;font-weight:bold;">PASS</span>'
        elif verdict == "REJECT":
            verdict_html = ' <span style="color:#e74c3c;font-weight:bold;">REJECT</span>'

        st.sidebar.markdown(
            f"""<div style="padding:6px 8px;margin-bottom:6px;border-left:3px solid {color};
            background:rgba(0,0,0,0.02);border-radius:0 4px 4px 0;font-size:13px;">
            {icon} <b>{agent}</b>{verdict_html}<br>
            <span style="color:#888;font-size:11px;">{seq}</span>
            </div>""",
            unsafe_allow_html=True,
        )

    # 전체 토론 기록 (접을 수 있는 영역)
    st.sidebar.markdown("---")
    st.sidebar.markdown("### Transcript Stream")

    # 에이전트별 통계
    pass_count = sum(1 for e in entries if e.get("verdict") == "PASS")
    reject_count = sum(1 for e in entries if e.get("verdict") == "REJECT")
    total_eval = pass_count + reject_count
    if total_eval > 0:
        st.sidebar.markdown(
            f"**Evaluations:** {total_eval} | "
            f'<span style="color:#2ecc71;">PASS {pass_count}</span> / '
            f'<span style="color:#e74c3c;">REJECT {reject_count}</span>',
            unsafe_allow_html=True,
        )

    # 최근 로그 항목 (최신순, 최대 20개)
    recent = list(reversed(entries[-20:]))
    for entry in recent:
        agent = entry["agent"]
        icon = AGENT_ICONS.get(agent, "🔹")
        color = AGENT_COLORS.get(agent, "#666")
        verdict = entry.get("verdict", "")
        seq = entry.get("sequence", "")
        msg = entry.get("message", "")

        # 메시지를 100자로 자르기
        short_msg = msg[:100] + "..." if len(msg) > 100 else msg

        verdict_tag = ""
        if verdict:
            vcolor = "#2ecc71" if verdict == "PASS" else "#e74c3c"
            verdict_tag = f' <span style="background:{vcolor};color:white;padding:1px 6px;border-radius:3px;font-size:10px;">{verdict}</span>'

        st.sidebar.markdown(
            f"""<div style="padding:4px 6px;margin-bottom:4px;border-left:2px solid {color};
            font-size:11px;background:rgba(0,0,0,0.01);border-radius:0 3px 3px 0;">
            {icon} <b>{agent}</b>{verdict_tag}
            <span style="color:#aaa;"> {seq}</span><br>
            <span style="color:#666;">{short_msg}</span>
            </div>""",
            unsafe_allow_html=True,
        )


# =========================================================================
#  UI 컴포넌트: Generated Library 테이블 (Module 00)
# =========================================================================

def render_raw_library_table():
    """Module 00 생성 서열 테이블 + 컬럼 설명을 렌더링한다."""
    df = load_raw_library()

    if df.empty:
        st.info("서열 라이브러리 대기 중... 00_peptide_gen.py 실행 후 표시됩니다.")
        return

    st.markdown("### Module 00 — Generated Peptide Library")

    # 요약 메트릭
    c1, c2, c3 = st.columns(3)
    c1.metric("Total Sequences", f"{len(df):,}")
    length_counts = df["Length"].value_counts().sort_index()
    c2.metric("Length Range", f"{df['Length'].min()}–{df['Length'].max()} AA")
    c3.metric("Unique Sequences", f"{df['Sequence'].nunique():,}")

    # 테이블 (검색 + 페이지네이션)
    col_s, col_n = st.columns([3, 1])
    with col_s:
        search_raw = st.text_input("Search sequence", "", key="search_raw",
                                   placeholder="e.g. FWKT, CPFW...")
    with col_n:
        rows_raw = st.selectbox("Rows", [20, 50, 100, 500], index=0, key="rows_raw")

    df_show = df.copy()
    if search_raw:
        df_show = df_show[df_show["Sequence"].str.contains(search_raw.upper(), na=False)]

    st.dataframe(
        df_show.head(rows_raw),
        use_container_width=True,
        hide_index=True,
        column_config={
            "Sequence": st.column_config.TextColumn("Sequence", width="large"),
            "Length": st.column_config.NumberColumn("Length (AA)", format="%d"),
        },
    )
    st.caption(f"Showing {min(rows_raw, len(df_show)):,} / {len(df_show):,} sequences"
               f"{' (filtered)' if search_raw else ''}")

    # 길이 분포 미니 차트
    with st.expander("Length Distribution", expanded=False):
        fig = px.histogram(df, x="Length", nbins=df["Length"].nunique(),
                           color_discrete_sequence=["#3498db"])
        fig.update_layout(
            xaxis_title="Sequence Length (AA)", yaxis_title="Count",
            template="plotly_white", height=280, margin=dict(t=20, b=40),
        )
        st.plotly_chart(fig, use_container_width=True)

    # 컬럼 설명 및 생성 기준
    st.markdown("""
    <div style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:8px;padding:14px 18px;margin-top:8px;font-size:13px;">
    <b>Column Descriptions & Generation Criteria</b>
    <table style="width:100%;margin-top:8px;border-collapse:collapse;">
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:4px 8px;font-weight:bold;width:140px;">Sequence</td>
      <td style="padding:4px 8px;">설계된 펩타이드 아미노산 서열. 모든 서열은 <b>FWKT</b>(Phe-Trp-Lys-Thr) 약리활성 모티프와
        <b>Cys-Cys</b> 이황화 결합 고리 구조를 필수 포함. Trp(W)은 SSTR2 소수성 포켓 삽입, Lys(K)는 Asp3.32 염교 형성 역할.</td>
    </tr>
    <tr>
      <td style="padding:4px 8px;font-weight:bold;">Length</td>
      <td style="padding:4px 8px;">서열 길이 (아미노산 수). <b>8–14 AA</b> 범위로 생성.
        Octreotide(8 AA)~Somatostatin-14(14 AA) 크기를 참조하여 SSTR2 결합 포켓에 적합한 길이로 제한.
        가변 영역 AA 풀: 17종 (표준 20종 − K, R, C). K는 FWKT 내부에서만 허용.</td>
    </tr>
    </table>
    </div>
    """, unsafe_allow_html=True)


# =========================================================================
#  UI 컴포넌트: PK-Filtered 테이블 (Module 01)
# =========================================================================

def render_filtered_table():
    """Module 01 PK 필터 통과 서열 테이블 + 컬럼 설명을 렌더링한다."""
    df = load_filtered()

    if df.empty:
        st.info("PK 필터 데이터 대기 중... 01_fast_pk_filter.py 실행 후 표시됩니다.")
        return

    st.markdown("### Module 01 — PK-Filtered Sequences")

    # Dotatate 기준값 (config.py에서 정의)
    DOTATATE_PI = 8.2212
    DOTATATE_HYDRO = 0.155

    # 요약 메트릭
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Passed Sequences", f"{len(df):,}")
    c2.metric("Avg Net Charge", f"{df['NetCharge'].mean():+.3f}")
    c3.metric("Avg pI", f"{df['pI'].mean():.2f}")
    c4.metric("Avg Hydrophobicity", f"{df['Hydrophobicity'].mean():.3f}")

    # 필터 컨트롤
    col_s, col_c, col_n = st.columns([3, 2, 1])
    with col_s:
        search_filt = st.text_input("Search sequence", "", key="search_filt",
                                    placeholder="e.g. CFWKT, CPFW...")
    with col_c:
        charge_filter = st.selectbox("Charge filter", [
            "All",
            f"Better than Dotatate (≤ {DOTATATE_CHARGE:+.3f})",
            "Negative charge only (< 0)",
        ], key="charge_filter")
    with col_n:
        rows_filt = st.selectbox("Rows", [20, 50, 100, 500], index=0, key="rows_filt")

    df_show = df.copy()
    if search_filt:
        df_show = df_show[df_show["Sequence"].str.contains(search_filt.upper(), na=False)]
    if "Better" in charge_filter:
        df_show = df_show[df_show["NetCharge"] <= DOTATATE_CHARGE]
    elif "Negative" in charge_filter:
        df_show = df_show[df_show["NetCharge"] < 0]

    # 컬러 하이라이트를 위해 스타일 적용
    def highlight_charge(val):
        if val <= 0:
            return "color: #2ecc71; font-weight: bold;"
        elif val <= DOTATATE_CHARGE:
            return "color: #27ae60;"
        else:
            return "color: #e67e22;"

    st.dataframe(
        df_show.head(rows_filt),
        use_container_width=True,
        hide_index=True,
        column_config={
            "Sequence": st.column_config.TextColumn("Sequence", width="large"),
            "Length": st.column_config.NumberColumn("Length", format="%d"),
            "NetCharge": st.column_config.NumberColumn("Net Charge", format="%+.4f"),
            "pI": st.column_config.NumberColumn("pI", format="%.4f"),
            "Hydrophobicity": st.column_config.NumberColumn("Hydrophobicity", format="%.4f"),
        },
    )
    st.caption(f"Showing {min(rows_filt, len(df_show)):,} / {len(df_show):,} sequences"
               f"{' (filtered)' if search_filt or 'All' not in charge_filter else ''}")

    # 분포 차트
    with st.expander("Property Distributions", expanded=False):
        fig_col1, fig_col2, fig_col3 = st.columns(3)
        with fig_col1:
            fig = px.histogram(df, x="NetCharge", nbins=40, color_discrete_sequence=["#2ecc71"])
            fig.add_vline(x=DOTATATE_CHARGE, line_dash="dash", line_color="#e74c3c",
                          annotation_text=f"Dotatate ({DOTATATE_CHARGE:+.3f})")
            fig.update_layout(xaxis_title="Net Charge", yaxis_title="Count",
                              template="plotly_white", height=250, margin=dict(t=20, b=40))
            st.plotly_chart(fig, use_container_width=True)
        with fig_col2:
            fig = px.histogram(df, x="pI", nbins=40, color_discrete_sequence=["#e67e22"])
            fig.add_vline(x=DOTATATE_PI, line_dash="dash", line_color="#e74c3c",
                          annotation_text=f"Dotatate ({DOTATATE_PI:.2f})")
            fig.update_layout(xaxis_title="pI", yaxis_title="Count",
                              template="plotly_white", height=250, margin=dict(t=20, b=40))
            st.plotly_chart(fig, use_container_width=True)
        with fig_col3:
            fig = px.histogram(df, x="Hydrophobicity", nbins=40, color_discrete_sequence=["#3498db"])
            fig.add_vline(x=DOTATATE_HYDRO, line_dash="dash", line_color="#e74c3c",
                          annotation_text=f"Dotatate ({DOTATATE_HYDRO:.3f})")
            fig.update_layout(xaxis_title="Hydrophobicity", yaxis_title="Count",
                              template="plotly_white", height=250, margin=dict(t=20, b=40))
            st.plotly_chart(fig, use_container_width=True)

    # 컬럼 설명 및 필터 기준
    st.markdown(f"""
    <div style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:8px;padding:14px 18px;margin-top:8px;font-size:13px;">
    <b>Column Descriptions & Filter Criteria</b>
    <table style="width:100%;margin-top:8px;border-collapse:collapse;">
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:5px 8px;font-weight:bold;width:160px;">Sequence</td>
      <td style="padding:5px 8px;">Module 00에서 생성되어 PK 필터를 통과한 펩타이드 서열.</td>
    </tr>
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:5px 8px;font-weight:bold;">Length</td>
      <td style="padding:5px 8px;">서열의 아미노산 수 (8–14 AA).</td>
    </tr>
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:5px 8px;font-weight:bold;">Net Charge</td>
      <td style="padding:5px 8px;">pH 7.4(혈액 pH)에서의 순전하. 양전하가 높을수록 신장 세뇨관의 메갈린(Megalin) 수용체에 흡착되어
        <b>신장 독성</b>이 증가한다.<br>
        <b>필터 기준: ≤ +1.0</b> (Dotatate 순전하 {DOTATATE_CHARGE:+.3f} 이하가 이상적).<br>
        <span style="color:#2ecc71;">● 음전하(−)</span> = 신장 안전성 우수 |
        <span style="color:#e67e22;">● 양전하(+)</span> = 주의 필요</td>
    </tr>
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:5px 8px;font-weight:bold;">pI</td>
      <td style="padding:5px 8px;">등전점 (Isoelectric Point). 순전하가 0이 되는 pH 값. pI &gt; 7이면 생리적 pH에서 양전하를 띄어
        신장 흡착 위험이 증가한다.<br>
        <b>필터 기준: ≤ {DOTATATE_PI:.2f} + 1.0 = ≤ 9.22</b> (Dotatate pI + 마진).</td>
    </tr>
    <tr>
      <td style="padding:5px 8px;font-weight:bold;">Hydrophobicity</td>
      <td style="padding:5px 8px;">Eisenberg 소수성 지수. 높을수록 지방 친화적(혈장 단백질 결합 ↑, 비특이적 조직 흡착 ↑),
        낮을수록 수용성이 좋다(혈중 반감기 ↑). 극단적 소수성은 신장 세뇨관 축적을 유발할 수 있다.<br>
        <b>필터 기준: Dotatate ({DOTATATE_HYDRO:.3f}) ± 0.5</b> (범위: {DOTATATE_HYDRO - 0.5:.3f} ~ {DOTATATE_HYDRO + 0.5:.3f}).</td>
    </tr>
    </table>
    </div>
    """, unsafe_allow_html=True)


# =========================================================================
#  UI 컴포넌트: 중앙 Candidate Ranking 테이블 (Module 02)
# =========================================================================

def render_candidate_table():
    """후보 랭킹 테이블을 3대 지표 PASS/FAIL과 함께 렌더링한다."""
    df = load_top100()

    if df.empty:
        st.info("후보 데이터 대기 중... 02_train_surrogate.py 실행 후 표시됩니다.")
        return

    st.markdown("## 📊 Candidate Ranking")

    # 필터 옵션
    col_f1, col_f2, col_f3 = st.columns([2, 2, 4])
    with col_f1:
        filter_mode = st.selectbox(
            "Filter",
            ["All", "All PASS", "Any FAIL", "Selectivity PASS", "Retention PASS", "Safety PASS"],
            key="filter_mode"
        )
    with col_f2:
        max_rows = st.selectbox("Show rows", [10, 25, 50, 100], index=0, key="max_rows")

    # 3대 지표 판정 컬럼 추가
    metrics_data = []
    for _, row in df.iterrows():
        m1, m2, m3 = judge_metrics(row)
        metrics_data.append({
            "Selectivity_PASS": m1,
            "Retention_PASS": m2,
            "Safety_PASS": m3,
            "All_PASS": m1 and m2 and m3,
        })
    df_metrics = pd.DataFrame(metrics_data)
    df_display = pd.concat([df.reset_index(drop=True), df_metrics], axis=1)

    # 필터링 적용
    if filter_mode == "All PASS":
        df_display = df_display[df_display["All_PASS"]]
    elif filter_mode == "Any FAIL":
        df_display = df_display[~df_display["All_PASS"]]
    elif filter_mode == "Selectivity PASS":
        df_display = df_display[df_display["Selectivity_PASS"]]
    elif filter_mode == "Retention PASS":
        df_display = df_display[df_display["Retention_PASS"]]
    elif filter_mode == "Safety PASS":
        df_display = df_display[df_display["Safety_PASS"]]

    # 요약 통계
    total = len(df)
    all_pass = df_metrics["All_PASS"].sum()
    sel_pass = df_metrics["Selectivity_PASS"].sum()
    ret_pass = df_metrics["Retention_PASS"].sum()
    saf_pass = df_metrics["Safety_PASS"].sum()

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("All 3 PASS", f"{all_pass}/{total}")
    c2.metric("① Selectivity ≥100x", f"{sel_pass}/{total}")
    c3.metric("② Retention (72h)", f"{ret_pass}/{total}")
    c4.metric("③ Safety (Charge)", f"{saf_pass}/{total}")

    # HTML 테이블 생성 (PASS/FAIL 뱃지 포함)
    html_rows = []
    for _, row in df_display.head(max_rows).iterrows():
        m1_badge = render_verdict_badge(row["Selectivity_PASS"])
        m2_badge = render_verdict_badge(row["Retention_PASS"])
        m3_badge = render_verdict_badge(row["Safety_PASS"])

        html_rows.append(
            f"<tr>"
            f"<td style='text-align:center;font-weight:bold;'>{int(row['Rank'])}</td>"
            f"<td><code style='background:#f0f0f0;padding:2px 6px;border-radius:3px;'>"
            f"{row['Sequence']}</code></td>"
            f"<td style='text-align:center;'>{row['Length']}</td>"
            f"<td style='text-align:center;font-weight:bold;color:#2c3e50;'>"
            f"{row['SSTR2_BindingEnergy']:.3f}</td>"
            f"<td style='text-align:center;'>{row['Selectivity']:.1f}x</td>"
            f"<td style='text-align:center;'>{row['CompositeScore']:.3f}</td>"
            f"<td style='text-align:center;'>{row['NetCharge']:+.4f}</td>"
            f"<td style='text-align:center;'>{m1_badge}</td>"
            f"<td style='text-align:center;'>{m2_badge}</td>"
            f"<td style='text-align:center;'>{m3_badge}</td>"
            f"</tr>"
        )

    table_html = f"""
    <div style="overflow-x:auto;">
    <table style="width:100%;border-collapse:collapse;font-size:13px;">
    <thead>
    <tr style="background:#2c3e50;color:white;">
        <th style="padding:8px 6px;">#</th>
        <th style="padding:8px 6px;">Sequence</th>
        <th style="padding:8px 6px;">Len</th>
        <th style="padding:8px 6px;">Binding Energy<br><span style="font-size:10px;">(kcal/mol)</span></th>
        <th style="padding:8px 6px;">Selectivity</th>
        <th style="padding:8px 6px;">Composite<br>Score</th>
        <th style="padding:8px 6px;">Net Charge</th>
        <th style="padding:8px 6px;">① Target<br>100x</th>
        <th style="padding:8px 6px;">② Retain<br>72h</th>
        <th style="padding:8px 6px;">③ Safety<br>Kidney</th>
    </tr>
    </thead>
    <tbody>
    {"".join(html_rows)}
    </tbody>
    </table>
    </div>
    """
    st.markdown(table_html, unsafe_allow_html=True)
    st.caption(f"Showing {min(max_rows, len(df_display))}/{len(df_display)} candidates (filtered from {total} total)")

    # 컬럼 설명 및 판정 기준
    st.markdown(f"""
    <div style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:8px;padding:14px 18px;margin-top:12px;font-size:13px;">
    <b>Column Descriptions & Success Criteria</b>
    <table style="width:100%;margin-top:8px;border-collapse:collapse;">
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:5px 8px;font-weight:bold;width:170px;">Binding Energy</td>
      <td style="padding:5px 8px;">XGBoost 대리 모델이 예측한 SSTR2 결합 자유 에너지 (kcal/mol).
        <b>음수값이 클수록 (더 음수일수록) 결합이 강하다</b>. Dotatate 기준값: {DOTATATE_DOCKING_SCORE} kcal/mol.</td>
    </tr>
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:5px 8px;font-weight:bold;">Selectivity</td>
      <td style="padding:5px 8px;">SSTR2 대비 비표적 수용체(SSTR1/3/4/5)에 대한 선택성 배수.
        Kd 비율(Kd<sub>non-target</sub> / Kd<sub>SSTR2</sub>)로 산출. <b>값이 높을수록 SSTR2에만 선택적으로 결합</b>한다.</td>
    </tr>
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:5px 8px;font-weight:bold;">Composite Score</td>
      <td style="padding:5px 8px;">결합력과 선택성을 종합한 순위 점수. <code>Binding Energy − log10(Selectivity)</code>로 산출.
        <b>더 음수일수록 우수</b>한 후보.</td>
    </tr>
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:5px 8px;font-weight:bold;">Net Charge</td>
      <td style="padding:5px 8px;">pH 7.4에서의 순전하. 양전하가 높으면 신장 독성 위험 증가. Dotatate: {DOTATATE_CHARGE:+.3f}.</td>
    </tr>
    </table>
    <div style="margin-top:10px;padding-top:8px;border-top:1px solid #dee2e6;">
    <b>3 Key Success Metrics (PASS/FAIL Criteria)</b>
    <table style="width:100%;margin-top:6px;border-collapse:collapse;">
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:4px 8px;font-weight:bold;width:170px;">① Target 100x</td>
      <td style="padding:4px 8px;"><b>Selectivity ≥ {SELECTIVITY_THRESHOLD:.0f}x</b> — SSTR2에 비표적 수용체 대비 100배 이상 선택적 결합.
        NET 종양의 SSTR2만 정밀 타격하여 off-target 부작용 최소화.</td>
    </tr>
    <tr style="border-bottom:1px solid #dee2e6;">
      <td style="padding:4px 8px;font-weight:bold;">② Retain 72h</td>
      <td style="padding:4px 8px;"><b>Binding Energy ≤ {MIN_RETENTION_ENERGY} kcal/mol</b> — 강한 결합 = 낮은 해리 속도(k<sub>off</sub>).
        k<sub>off</sub> &lt; 1/(72h) = 3.86×10⁻⁶ /s 수준의 체류력으로 Lu-177 방사선 조사 시간 확보.</td>
    </tr>
    <tr>
      <td style="padding:4px 8px;font-weight:bold;">③ Safety Kidney</td>
      <td style="padding:4px 8px;"><b>Net Charge ≤ {DOTATATE_CHARGE:+.3f}</b> (Dotatate 이하) — 양전하 감소 시 메갈린(Megalin) 매개
        신장 세뇨관 재흡수가 감소하여 신장 방사선 피폭 20% 이상 개선 기대.</td>
    </tr>
    </table>
    </div>
    </div>
    """, unsafe_allow_html=True)


# =========================================================================
#  UI 컴포넌트: Dotatate 비교 차트
# =========================================================================

def render_comparison_charts():
    """후보물질 vs Dotatate 비교 시각화를 렌더링한다."""
    df = load_top100()
    if df.empty:
        return

    st.markdown("## 📈 Performance Analysis")

    tab1, tab2, tab3 = st.tabs([
        "Binding Energy Distribution",
        "Selectivity vs Binding",
        "Net Charge Comparison",
    ])

    # --- Tab 1: Binding Energy Distribution ---
    with tab1:
        fig = go.Figure()
        fig.add_trace(go.Histogram(
            x=df["SSTR2_BindingEnergy"],
            nbinsx=30,
            name="Top Candidates",
            marker_color="#3498db",
            opacity=0.8,
        ))
        fig.add_vline(
            x=DOTATATE_DOCKING_SCORE,
            line_dash="dash", line_color="#e74c3c", line_width=2,
            annotation_text=f"Dotatate ({DOTATATE_DOCKING_SCORE})",
            annotation_position="top right",
        )
        fig.add_vline(
            x=MIN_RETENTION_ENERGY,
            line_dash="dot", line_color="#f39c12", line_width=2,
            annotation_text=f"72h Retention ({MIN_RETENTION_ENERGY})",
            annotation_position="top left",
        )
        fig.update_layout(
            title="SSTR2 Binding Energy Distribution (kcal/mol)",
            xaxis_title="Binding Energy (kcal/mol)",
            yaxis_title="Count",
            template="plotly_white",
            height=400,
        )
        st.plotly_chart(fig, use_container_width=True)

    # --- Tab 2: Selectivity vs Binding Scatter ---
    with tab2:
        fig = px.scatter(
            df, x="SSTR2_BindingEnergy", y="Selectivity",
            color="CompositeScore",
            size=np.clip(df["Selectivity"] / df["Selectivity"].max() * 15, 3, 15).tolist(),
            hover_data=["Rank", "Sequence", "NetCharge"],
            color_continuous_scale="Viridis",
        )
        # 성공 기준 영역 표시
        fig.add_hline(
            y=SELECTIVITY_THRESHOLD, line_dash="dash", line_color="#e74c3c",
            annotation_text="100x Selectivity Threshold",
        )
        fig.add_vline(
            x=MIN_RETENTION_ENERGY, line_dash="dash", line_color="#f39c12",
            annotation_text="72h Retention Threshold",
        )
        # "Winning zone" 사각형
        fig.add_shape(
            type="rect",
            x0=df["SSTR2_BindingEnergy"].min() - 0.5, x1=MIN_RETENTION_ENERGY,
            y0=SELECTIVITY_THRESHOLD, y1=df["Selectivity"].max() * 1.1,
            fillcolor="rgba(46,204,113,0.08)", line=dict(color="#2ecc71", width=2, dash="dot"),
        )
        fig.add_annotation(
            x=(df["SSTR2_BindingEnergy"].min() + MIN_RETENTION_ENERGY) / 2,
            y=df["Selectivity"].max() * 0.95,
            text="<b>WINNING ZONE</b>",
            font=dict(color="#2ecc71", size=14),
            showarrow=False,
        )
        fig.update_layout(
            title="Selectivity vs. Binding Affinity",
            xaxis_title="SSTR2 Binding Energy (kcal/mol)",
            yaxis_title="Selectivity (fold)",
            template="plotly_white",
            height=450,
        )
        st.plotly_chart(fig, use_container_width=True)

    # --- Tab 3: Net Charge Comparison ---
    with tab3:
        # Top 10 후보의 전하 vs Dotatate
        df_top10 = df.head(10).copy()
        fig = go.Figure()
        # Dotatate 기준선
        fig.add_trace(go.Bar(
            x=["Dotatate (ref)"],
            y=[DOTATATE_CHARGE],
            name="Dotatate",
            marker_color="#e74c3c",
            width=0.4,
        ))
        fig.add_trace(go.Bar(
            x=[f"#{int(r['Rank'])} {r['Sequence'][:6]}.." for _, r in df_top10.iterrows()],
            y=df_top10["NetCharge"].tolist(),
            name="Candidates",
            marker_color=["#2ecc71" if c <= DOTATATE_CHARGE else "#f39c12"
                          for c in df_top10["NetCharge"]],
            width=0.4,
        ))
        fig.add_hline(
            y=DOTATATE_CHARGE, line_dash="dash", line_color="#e74c3c",
            annotation_text=f"Dotatate Charge ({DOTATATE_CHARGE:+.3f})",
        )
        fig.update_layout(
            title="Net Charge Comparison: Top 10 Candidates vs. Dotatate",
            xaxis_title="Peptide",
            yaxis_title="Net Charge (pH 7.4)",
            template="plotly_white",
            height=400,
            showlegend=False,
        )
        st.plotly_chart(fig, use_container_width=True)


# =========================================================================
#  UI 컴포넌트: 최종 보고서 뷰어
# =========================================================================

def render_report_viewer():
    """최종 보고서 마크다운을 렌더링한다."""
    report = load_report()
    if report:
        st.markdown("## 📄 Final Winning Peptides Report")
        with st.expander("View Full Report", expanded=False):
            st.markdown(report, unsafe_allow_html=True)
    else:
        st.info("최종 보고서 대기 중... 03_ai_boardroom.py 완료 후 표시됩니다.")


# =========================================================================
#  페이지 설정 및 메인 레이아웃
# =========================================================================

def main():
    st.set_page_config(
        page_title="SSTR2 AI Scientist - Mission Control",
        page_icon="🧪",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    # 커스텀 CSS
    st.markdown("""
    <style>
    /* 헤더 스타일 */
    .main-header {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        padding: 20px 30px;
        border-radius: 10px;
        margin-bottom: 20px;
        color: white;
    }
    .main-header h1 {
        margin: 0;
        font-size: 28px;
        font-weight: 700;
    }
    .main-header p {
        margin: 5px 0 0 0;
        opacity: 0.8;
        font-size: 14px;
    }
    /* 테이블 행 호버 */
    table tbody tr:hover { background: rgba(52, 152, 219, 0.08) !important; }
    table tbody tr:nth-child(even) { background: rgba(0,0,0,0.02); }
    table th, table td { padding: 8px 6px !important; border-bottom: 1px solid #eee; }
    /* metric 카드 */
    [data-testid="stMetric"] {
        background: rgba(0,0,0,0.02);
        border-radius: 8px;
        padding: 10px;
        border: 1px solid #eee;
    }
    /* 사이드바 */
    [data-testid="stSidebar"] { background: #fafbfc; }
    </style>
    """, unsafe_allow_html=True)

    # ---- 헤더 ----
    st.markdown(
        """<div class="main-header">
        <h1>🧪 SSTR2 AI Scientist — Mission Control</h1>
        <p>Real-time Pipeline Monitor | 3 Key Metrics: 100x Selectivity · 72h Retention · 20% Kidney Safety</p>
        </div>""",
        unsafe_allow_html=True,
    )

    # ---- 사이드바: Agent Monitor ----
    render_agent_monitor()

    # 사이드바 하단: 자동 새로고침
    st.sidebar.markdown("---")
    auto_refresh = st.sidebar.checkbox("Auto-refresh (5s)", value=False)
    if st.sidebar.button("🔄 Refresh Now"):
        st.cache_data.clear()
        st.rerun()

    # ---- 상단: Pipeline Progress ----
    render_pipeline_progress()
    st.markdown("---")

    # ---- 중앙: 3개 데이터셋 탭 ----
    st.markdown("## 📊 Data Explorer")
    tab_top100, tab_filtered, tab_raw = st.tabs([
        "🏆 Top-100 Candidates (Module 02)",
        "🔬 PK-Filtered Sequences (Module 01)",
        "🧬 Generated Library (Module 00)",
    ])
    with tab_top100:
        render_candidate_table()
    with tab_filtered:
        render_filtered_table()
    with tab_raw:
        render_raw_library_table()

    st.markdown("---")

    # ---- 하단: 분석 차트 ----
    render_comparison_charts()
    st.markdown("---")

    # ---- 최종 보고서 ----
    render_report_viewer()

    # ---- 자동 새로고침 ----
    if auto_refresh:
        time.sleep(5)
        st.cache_data.clear()
        st.rerun()


if __name__ == "__main__":
    main()
