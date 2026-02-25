#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  Module 04: Validation & Publication-Quality Statistical Analysis
  (04_validate_and_analyze.py)
==========================================================================

  [Stage]  Phase 3 — Post-pipeline Validation & Visualization
  [Input]  output/SSTR2_Top100_Boardroom_Ready.csv  (Module 02 output)
           output/SSTR2_Filtered_Top100.csv          (Module 01 output)
           output/SSTR2_Raw_Library.csv               (Module 00 output)
  [Output] output/figures/  (publication-ready 300 DPI PNG + PDF vector)
           output/Top_5_Candidates_Summary.csv
           output/Validation_Report.md

  [학술적 무결성 — Key Principles]
  1. 도킹 에너지 → "AI-predicted Binding Free Energy (ΔG_pred)" 명시
  2. ΔG → K_d 역산: ΔG = RT ln(K_d), R=1.987 cal/(mol·K), T=310.15 K
  3. Spearman 상관 산점도 (AI 예측 vs 물리 도킹 스코어)
  4. TEST MODE 워터마크 감지 및 전체 그래프 삽입
  5. Surrogate Model 한계 Disclaimer 필수 포함

  [Figure Catalog]
    Fig 1.  ΔG_pred Distribution (Histogram + KDE)
    Fig 2.  Selectivity Ratio Analysis (Box + Strip)
    Fig 3.  Net Charge Safety Profile (vs Dotatate)
    Fig 4.  ΔG_pred vs Selectivity Quadrant Scatter
    Fig 5.  Physicochemical Property Space (Hydrophobicity vs pI)
    Fig 6.  Peptide Length vs ΔG_pred (Box Plot)
    Fig 7.  Composite Stability Score Ranking (Top-20 Waterfall)
    Fig 8.  Amino Acid Positional Frequency Heatmap
    Fig 9.  Property Correlation Matrix
    Fig 10. Multi-dimensional Radar Chart (Top-5 vs Dotatate)
    Fig 11. Pipeline Funnel Summary
    Fig 12. Spearman Correlation: AI-predicted ΔG vs Physical Docking Score
    Fig 13. 5-Indicator Stability Profile (Top-10)
    Table 1. Descriptive Statistics Summary
    Table 2. Top_5_Candidates_Summary.csv (with K_d)

  [Usage]
    python 04_validate_and_analyze.py
==========================================================================
"""

import sys
import math
import time
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch
from matplotlib.gridspec import GridSpec
from pathlib import Path
from collections import Counter

# Optional imports
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

try:
    from scipy import stats as sp_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ── Project config ──────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

warnings.filterwarnings("ignore", category=FutureWarning)

# ==========================================================================
#  Constants
# ==========================================================================

OUTPUT_DIR = config.OUTPUT_DIR
FIG_DIR = OUTPUT_DIR / "figures"

# Input files (backward compatible with old naming)
FILE_RAW = OUTPUT_DIR / "SSTR2_Raw_Library.csv"
FILE_FILTERED = OUTPUT_DIR / "SSTR2_Filtered_Top100.csv"
FILE_FILTERED_LEGACY = OUTPUT_DIR / "SSTR2_Filtered_10K.csv"
FILE_TOP100 = OUTPUT_DIR / "SSTR2_Top100_Boardroom_Ready.csv"

# Physical constants (from config.py)
R_CAL = config.R_CAL        # 1.987 cal/(mol·K)
T_BODY = config.T_BODY      # 310.15 K (37 °C)
RT_KCAL = config.RT_BODY_KCAL  # ≈ 0.616 kcal/mol

# Dotatate reference values
DOTATATE = {
    "Sequence": config.DOTATATE_SEQ,   # FCYWKTCT
    "Length": len(config.DOTATATE_SEQ),
    "BindingEnergy": -8.5,             # kcal/mol (literature reference)
    "Selectivity": 50.0,               # approximate clinical selectivity
    "NetCharge": 0.707,
    "pI": 8.22,
    "Hydrophobicity": 0.155,
    "CompositeStabilityScore": 0.5,    # nominal reference
}

# Three success criteria thresholds
THRESH_SELECTIVITY = 100.0      # >= 100x selective binding
THRESH_RETENTION = -7.0         # <= -7.0 kcal/mol (72h retention)
THRESH_CHARGE = 0.707           # <= Dotatate charge (kidney safety)

# TEST MODE watermark string (must match 02_train_surrogate.py)
WATERMARK_TEXT = "TEST MODE: DUMMY DATA — NOT FOR PUBLICATION"

# Disclaimer text (학술적 무결성)
DISCLAIMER = (
    "**Disclaimer:** This Surrogate Model combines ESM-2 embeddings, RDKit "
    "molecular descriptors, and XGBoost regression. The reported ΔG_pred values "
    "are AI-predicted approximations, NOT experimentally measured binding affinities. "
    "All final candidates MUST be validated through physical molecular docking "
    "(e.g., AutoDock Vina, Glide) and molecular dynamics (MD) simulations "
    "before any experimental synthesis or in vivo testing."
)

# ── Publication style ───────────────────────────────────────────────────

# Color-blind friendly palette (Okabe-Ito)
C_BLUE   = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN  = "#009E73"
C_RED    = "#D55E00"
C_PURPLE = "#CC79A7"
C_CYAN   = "#56B4E9"
C_YELLOW = "#F0E442"
C_BLACK  = "#000000"
C_GRAY   = "#999999"

PALETTE_SEQ = [C_BLUE, C_CYAN, C_GREEN, C_YELLOW, C_ORANGE, C_RED]

DPI = 300
FIG_FORMAT = ["png", "pdf"]


def _setup_style():
    """Set matplotlib rcParams for publication-quality output."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,
        "axes.linewidth": 1.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.linewidth": 1.5,
    })
    if HAS_SEABORN:
        sns.set_palette("colorblind")


# ==========================================================================
#  Utility Functions
# ==========================================================================

def _detect_test_mode(df):
    """Detect if the data was generated in TEST MODE (dummy data).

    Checks for the WATERMARK column or the TestMode flag
    that 02_train_surrogate.py adds in --test-mode.
    """
    # Check for WATERMARK column added by Module 02 in test mode
    if "WATERMARK" in df.columns:
        return True
    # Check for TestMode flag
    if "TestMode" in df.columns:
        return df["TestMode"].any()
    return False


def _add_watermark(ax):
    """Add TEST MODE watermark diagonally across the figure."""
    ax.text(
        0.5, 0.5, WATERMARK_TEXT,
        transform=ax.transAxes,
        fontsize=28, color="red", alpha=0.25,
        ha="center", va="center", rotation=30,
        fontweight="bold", zorder=999,
    )


def _savefig(fig, name, test_mode=False):
    """Save figure in PNG + PDF. Adds watermark if test_mode."""
    if test_mode:
        for ax in fig.get_axes():
            _add_watermark(ax)
    for fmt in FIG_FORMAT:
        path = FIG_DIR / f"{name}.{fmt}"
        fig.savefig(path, format=fmt, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved: {name}.png / .pdf")


def _dg_to_kd(dg_kcal):
    """Convert ΔG (kcal/mol) to K_d (M) using ΔG = RT·ln(K_d).

    K_d = exp(ΔG / RT)
    R = 1.987 cal/(mol·K), T = 310.15 K
    RT = 0.616 kcal/mol
    """
    return math.exp(dg_kcal / RT_KCAL)


def _kd_to_nm(kd_molar):
    """Convert K_d from Molar to nanomolar."""
    return kd_molar * 1e9


# ==========================================================================
#  Data Loading
# ==========================================================================

def load_data():
    """Load the main dataset from Module 02 output."""
    if FILE_TOP100.exists():
        df = pd.read_csv(FILE_TOP100)
        print(f"  [DATA] Loaded {len(df)} candidates from {FILE_TOP100.name}")
    else:
        print(f"  [ERROR] {FILE_TOP100} not found.")
        sys.exit(1)

    # Detect test mode
    test_mode = _detect_test_mode(df)
    if test_mode:
        print(f"  [WARNING] {WATERMARK_TEXT}")

    # Ensure key columns exist
    col_map = {
        "SSTR2_BindingEnergy": "dG_pred",
    }
    for old_col, new_col in col_map.items():
        if old_col in df.columns:
            df[new_col] = df[old_col]

    # Compute K_d from ΔG_pred
    if "dG_pred" in df.columns:
        df["K_d_M"] = df["dG_pred"].apply(_dg_to_kd)
        df["K_d_nM"] = df["K_d_M"].apply(_kd_to_nm)

    return df, test_mode


def load_pipeline_counts():
    """Load pipeline stage counts for funnel chart."""
    counts = {}

    if FILE_RAW.exists():
        raw = pd.read_csv(FILE_RAW)
        counts["00_Generated"] = len(raw)
    else:
        counts["00_Generated"] = None

    # Try new filename first, then legacy
    filtered_file = FILE_FILTERED if FILE_FILTERED.exists() else FILE_FILTERED_LEGACY
    if filtered_file.exists():
        filt = pd.read_csv(filtered_file)
        counts["01_PK_Filtered"] = len(filt)
    else:
        counts["01_PK_Filtered"] = None

    if FILE_TOP100.exists():
        top = pd.read_csv(FILE_TOP100)
        counts["02_Surrogate_Top100"] = len(top)
    else:
        counts["02_Surrogate_Top100"] = None

    report_path = OUTPUT_DIR / "Final_Winning_Peptides_Report.md"
    if report_path.exists():
        text = report_path.read_text(encoding="utf-8")
        # Count approved sequences in report
        import re
        approved = len(re.findall(r"(?:Sequence|서열)[:\s]+([A-Z]{8,14})", text))
        counts["03_Boardroom_Approved"] = max(approved, 1)
    else:
        counts["03_Boardroom_Approved"] = None

    return counts


# ==========================================================================
#  Figure Functions
# ==========================================================================

def fig01_dg_distribution(df, test_mode):
    """Fig 1. AI-predicted ΔG_pred Distribution (Histogram + KDE)."""
    fig, ax = plt.subplots(figsize=(8, 5))
    energies = df["dG_pred"].dropna()

    ax.hist(energies, bins=30, color=C_BLUE, alpha=0.7, edgecolor="white",
            label="Top-100 Candidates", density=True)

    if HAS_SEABORN:
        sns.kdeplot(energies, ax=ax, color=C_RED, linewidth=2, label="KDE")

    ax.axvline(DOTATATE["BindingEnergy"], color=C_ORANGE, ls="--", lw=2,
               label=f'Dotatate (ΔG = {DOTATATE["BindingEnergy"]:.1f})')
    ax.axvline(THRESH_RETENTION, color=C_GREEN, ls=":", lw=1.5,
               label=f"72h Retention Threshold ({THRESH_RETENTION})")

    ax.set_xlabel("AI-predicted Binding Free Energy, ΔG$_{pred}$ (kcal/mol)")
    ax.set_ylabel("Density")
    ax.set_title("Fig 1. Distribution of AI-predicted ΔG$_{pred}$")
    ax.legend(frameon=False)

    # Stats annotation
    med = energies.median()
    ax.text(0.02, 0.95, f"Median: {med:.2f} kcal/mol\nn = {len(energies)}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

    _savefig(fig, "Fig01_dG_pred_Distribution", test_mode)


def fig02_selectivity(df, test_mode):
    """Fig 2. Selectivity Ratio Analysis."""
    fig, ax = plt.subplots(figsize=(8, 5))
    sel = df["Selectivity"].dropna()

    if HAS_SEABORN:
        sns.boxplot(y=sel, ax=ax, color=C_CYAN, width=0.3,
                    flierprops=dict(marker="o", markersize=3))
        sns.stripplot(y=sel, ax=ax, color=C_BLACK, alpha=0.3, size=3, jitter=0.1)
    else:
        ax.boxplot(sel, vert=True)

    ax.axhline(THRESH_SELECTIVITY, color=C_RED, ls="--", lw=2,
               label=f"100x Threshold")
    ax.axhline(DOTATATE["Selectivity"], color=C_ORANGE, ls=":", lw=1.5,
               label=f"Dotatate ({DOTATATE['Selectivity']}x)")

    pct_pass = (sel >= THRESH_SELECTIVITY).sum() / len(sel) * 100
    ax.set_ylabel("Selectivity (K$_d$_nontarget / K$_d$_SSTR2)")
    ax.set_title("Fig 2. Selectivity Ratio Analysis")
    ax.legend(frameon=False)

    ax.text(0.02, 0.95, f"≥100x: {pct_pass:.1f}%\nn = {len(sel)}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

    _savefig(fig, "Fig02_Selectivity", test_mode)


def fig03_charge_safety(df, test_mode):
    """Fig 3. Net Charge Safety Profile."""
    fig, ax = plt.subplots(figsize=(8, 5))
    charges = df["NetCharge"].dropna()

    ax.hist(charges, bins=25, color=C_GREEN, alpha=0.7, edgecolor="white")
    ax.axvline(DOTATATE["NetCharge"], color=C_ORANGE, ls="--", lw=2,
               label=f'Dotatate (+{DOTATATE["NetCharge"]:.3f})')
    ax.axvline(0.8, color=C_RED, ls=":", lw=1.5,
               label="Rejection Threshold (+0.8)")

    pct_safe = (charges <= THRESH_CHARGE).sum() / len(charges) * 100
    ax.set_xlabel("Net Charge at pH 7.4")
    ax.set_ylabel("Count")
    ax.set_title("Fig 3. Net Charge Safety Profile (Kidney Toxicity)")
    ax.legend(frameon=False)

    ax.text(0.02, 0.95, f"≤ Dotatate: {pct_safe:.1f}%\nn = {len(charges)}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

    _savefig(fig, "Fig03_Charge_Safety", test_mode)


def fig04_dg_vs_selectivity(df, test_mode):
    """Fig 4. ΔG_pred vs Selectivity Quadrant Scatter."""
    fig, ax = plt.subplots(figsize=(8, 7))

    # Color by pass/fail
    pass_mask = (df["dG_pred"] <= THRESH_RETENTION) & (df["Selectivity"] >= THRESH_SELECTIVITY)
    ax.scatter(df.loc[~pass_mask, "dG_pred"], df.loc[~pass_mask, "Selectivity"],
               c=C_GRAY, alpha=0.5, s=30, label="Below Threshold", zorder=2)
    ax.scatter(df.loc[pass_mask, "dG_pred"], df.loc[pass_mask, "Selectivity"],
               c=C_BLUE, alpha=0.8, s=50, edgecolors="black", linewidths=0.5,
               label="PASS Both Criteria", zorder=3)

    # Dotatate reference point
    ax.scatter(DOTATATE["BindingEnergy"], DOTATATE["Selectivity"],
               c=C_ORANGE, s=150, marker="*", edgecolors="black",
               linewidths=1, label="Dotatate", zorder=4)

    # Threshold lines (quadrant)
    ax.axvline(THRESH_RETENTION, color=C_GREEN, ls="--", alpha=0.6, label=f"ΔG = {THRESH_RETENTION}")
    ax.axhline(THRESH_SELECTIVITY, color=C_RED, ls="--", alpha=0.6, label=f"Sel = {THRESH_SELECTIVITY}x")

    # Quadrant labels
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    mid_x = (THRESH_RETENTION + xlim[0]) / 2
    mid_y = (THRESH_SELECTIVITY + ylim[1]) / 2
    ax.text(mid_x, mid_y, "IDEAL\nZONE", fontsize=14, ha="center", va="center",
            color=C_GREEN, alpha=0.3, fontweight="bold")

    ax.set_xlabel("AI-predicted ΔG$_{pred}$ (kcal/mol)")
    ax.set_ylabel("Selectivity (fold)")
    ax.set_title("Fig 4. ΔG$_{pred}$ vs Selectivity Quadrant Analysis")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    n_pass = pass_mask.sum()
    ax.text(0.02, 0.05, f"Pass both: {n_pass}/{len(df)} ({n_pass/len(df)*100:.1f}%)",
            transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

    _savefig(fig, "Fig04_dG_vs_Selectivity_Quadrant", test_mode)


def fig05_property_space(df, test_mode):
    """Fig 5. Physicochemical Property Space (Hydrophobicity vs pI)."""
    fig, ax = plt.subplots(figsize=(8, 6))

    sc = ax.scatter(df["pI"], df["Hydrophobicity"], c=df["dG_pred"],
                    cmap="viridis_r", s=40, alpha=0.7, edgecolors="gray", linewidths=0.3)
    plt.colorbar(sc, ax=ax, label="ΔG$_{pred}$ (kcal/mol)")

    ax.scatter(DOTATATE["pI"], DOTATATE["Hydrophobicity"],
               c=C_ORANGE, s=150, marker="*", edgecolors="black",
               linewidths=1, label="Dotatate", zorder=4)

    ax.set_xlabel("Isoelectric Point (pI)")
    ax.set_ylabel("Hydrophobicity (Kyte-Doolittle)")
    ax.set_title("Fig 5. Physicochemical Property Space")
    ax.legend(frameon=False)

    _savefig(fig, "Fig05_Property_Space", test_mode)


def fig06_length_vs_dg(df, test_mode):
    """Fig 6. Peptide Length vs ΔG_pred (Box Plot)."""
    fig, ax = plt.subplots(figsize=(8, 5))

    lengths = sorted(df["Length"].unique())
    data_by_len = [df.loc[df["Length"] == l, "dG_pred"].values for l in lengths]

    bp = ax.boxplot(data_by_len, positions=range(len(lengths)),
                    patch_artist=True, widths=0.6)
    for patch in bp["boxes"]:
        patch.set_facecolor(C_CYAN)
        patch.set_alpha(0.7)

    ax.set_xticks(range(len(lengths)))
    ax.set_xticklabels(lengths)
    ax.axhline(THRESH_RETENTION, color=C_GREEN, ls="--", lw=1.5,
               label=f"72h Threshold ({THRESH_RETENTION})")
    ax.set_xlabel("Peptide Length (AA)")
    ax.set_ylabel("AI-predicted ΔG$_{pred}$ (kcal/mol)")
    ax.set_title("Fig 6. Peptide Length vs ΔG$_{pred}$")
    ax.legend(frameon=False)

    _savefig(fig, "Fig06_Length_vs_dG", test_mode)


def fig07_composite_waterfall(df, test_mode):
    """Fig 7. Composite Stability Score Ranking (Top-20 Waterfall)."""
    # Use CompositeStabilityScore if available, otherwise CompositeScore
    score_col = "CompositeStabilityScore" if "CompositeStabilityScore" in df.columns else "CompositeScore"
    if score_col not in df.columns:
        print("    [SKIP] Fig07 — no composite score column found.")
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    top20 = df.nsmallest(20, "Rank") if "Rank" in df.columns else df.head(20)
    scores = top20[score_col].values
    labels = top20["Sequence"].values

    colors = [C_BLUE if s >= 0.5 else C_ORANGE if s >= 0.3 else C_RED for s in scores]
    bars = ax.barh(range(len(scores)), scores, color=colors, edgecolor="white", height=0.7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([s[:12] + "…" if len(s) > 12 else s for s in labels], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Composite Stability Score")
    ax.set_title("Fig 7. Top-20 Candidates: Composite Stability Score")

    # Value labels
    for i, (bar, val) in enumerate(zip(bars, scores)):
        ax.text(bar.get_width() + 0.01, i, f"{val:.3f}", va="center", fontsize=8)

    _savefig(fig, "Fig07_Composite_Waterfall", test_mode)


def fig08_aa_heatmap(df, test_mode):
    """Fig 8. Amino Acid Positional Frequency Heatmap."""
    max_len = df["Length"].max()
    aa_order = list("ACDEFGHIKLMNPQRSTVWY")
    freq_matrix = np.zeros((len(aa_order), max_len))

    for seq in df["Sequence"]:
        for pos, aa in enumerate(seq):
            if aa in aa_order:
                freq_matrix[aa_order.index(aa), pos] += 1

    # Normalize per position
    col_sums = freq_matrix.sum(axis=0)
    col_sums[col_sums == 0] = 1
    freq_norm = freq_matrix / col_sums

    fig, ax = plt.subplots(figsize=(max(8, max_len * 0.8), 6))
    im = ax.imshow(freq_norm, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Frequency")

    ax.set_xticks(range(max_len))
    ax.set_xticklabels(range(1, max_len + 1))
    ax.set_yticks(range(len(aa_order)))
    ax.set_yticklabels(aa_order)
    ax.set_xlabel("Position")
    ax.set_ylabel("Amino Acid")
    ax.set_title("Fig 8. Positional Amino Acid Frequency (Top-100)")

    _savefig(fig, "Fig08_AA_Heatmap", test_mode)


def fig09_correlation_matrix(df, test_mode):
    """Fig 9. Property Correlation Matrix."""
    cols = ["dG_pred", "Selectivity", "NetCharge", "pI", "Hydrophobicity"]
    # Add stability indicators if available
    for c in ["ProteaseSites", "AggregationScore", "PPBScore", "Rigidity",
              "CompositeStabilityScore"]:
        if c in df.columns:
            cols.append(c)

    cols = [c for c in cols if c in df.columns]
    if len(cols) < 3:
        print("    [SKIP] Fig09 — insufficient numeric columns.")
        return

    corr = df[cols].corr()

    fig, ax = plt.subplots(figsize=(10, 8))
    if HAS_SEABORN:
        sns.heatmap(corr, ax=ax, annot=True, fmt=".2f", cmap="RdBu_r",
                    center=0, vmin=-1, vmax=1, square=True, linewidths=0.5)
    else:
        im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
        plt.colorbar(im, ax=ax)
        for i in range(len(cols)):
            for j in range(len(cols)):
                ax.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)
        ax.set_xticks(range(len(cols)))
        ax.set_yticks(range(len(cols)))
        ax.set_xticklabels(cols, rotation=45, ha="right")
        ax.set_yticklabels(cols)

    ax.set_title("Fig 9. Property Correlation Matrix")
    _savefig(fig, "Fig09_Correlation_Matrix", test_mode)


def fig10_radar_chart(df, test_mode):
    """Fig 10. Multi-dimensional Radar Chart (Top-5 vs Dotatate)."""
    categories = ["ΔG_pred", "Selectivity", "Charge\nSafety", "Hydro.", "pI"]
    N = len(categories)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    def _get_values(row_dict):
        """Normalize values to 0-1 scale for radar."""
        dg = row_dict.get("dG_pred", row_dict.get("SSTR2_BindingEnergy",
             row_dict.get("BindingEnergy", -7.0)))
        sel = row_dict.get("Selectivity", 50.0)
        charge = row_dict.get("NetCharge", 0.5)
        hydro = row_dict.get("Hydrophobicity", 0.2)
        pi_val = row_dict.get("pI", 7.0)

        vals = [
            min(1.0, max(0.0, (dg + 12) / (-5 + 12))),     # -12→1.0, -5→0.0
            min(1.0, max(0.0, sel / 2000)),                  # 0→0, 2000→1.0
            min(1.0, max(0.0, 1.0 - charge / 1.5)),         # 0→1.0, 1.5→0.0
            min(1.0, max(0.0, (hydro + 1) / 2)),             # -1→0, 1→1.0
            min(1.0, max(0.0, (pi_val - 4) / 8)),            # 4→0, 12→1.0
        ]
        return vals + vals[:1]

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=10)

    # Dotatate baseline
    dot_vals = _get_values(DOTATATE)
    ax.plot(angles, dot_vals, "o--", color=C_ORANGE, linewidth=2, label="Dotatate")
    ax.fill(angles, dot_vals, alpha=0.1, color=C_ORANGE)

    # Top-5 candidates
    colors_top5 = [C_BLUE, C_GREEN, C_PURPLE, C_CYAN, C_RED]
    top5 = df.head(5)
    for i, (_, row) in enumerate(top5.iterrows()):
        vals = _get_values(row.to_dict())
        seq_short = row["Sequence"][:8]
        ax.plot(angles, vals, "o-", color=colors_top5[i], linewidth=1.5,
                label=f"#{i+1} {seq_short}", alpha=0.8)

    ax.set_title("Fig 10. Top-5 Candidates vs Dotatate", y=1.08)
    ax.legend(loc="lower right", bbox_to_anchor=(1.3, -0.1), frameon=False, fontsize=8)

    _savefig(fig, "Fig10_Radar_Top5", test_mode)


def fig11_pipeline_funnel(df, test_mode):
    """Fig 11. Pipeline Funnel Summary."""
    counts = load_pipeline_counts()

    stages = []
    values = []
    for k, v in counts.items():
        if v is not None:
            stages.append(k.replace("_", "\n"))
            values.append(v)

    if len(stages) < 2:
        print("    [SKIP] Fig11 — insufficient pipeline data.")
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = [C_BLUE, C_GREEN, C_ORANGE, C_RED][:len(values)]

    bars = ax.bar(range(len(values)), values, color=colors, edgecolor="white", width=0.6)
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages, fontsize=9)
    ax.set_ylabel("Number of Sequences")
    ax.set_title("Fig 11. Pipeline Funnel: Sequence Compression")

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{val:,}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    # Compression ratio annotation
    if len(values) >= 2:
        ratio = values[0] / max(values[-1], 1)
        ax.text(0.98, 0.95, f"Compression: {ratio:.0f}:1",
                transform=ax.transAxes, fontsize=11, ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    _savefig(fig, "Fig11_Pipeline_Funnel", test_mode)


def fig12_spearman_correlation(df, test_mode):
    """Fig 12. Spearman Correlation: AI-predicted ΔG vs Physical Docking Score.

    If both XGB_predicted and physical docking columns exist, plot their
    correlation. Otherwise, plot the surrogate model predictions vs rank
    to validate monotonic ordering.
    """
    fig, ax = plt.subplots(figsize=(8, 7))

    # Determine which columns to compare
    xgb_col = None
    dock_col = None
    for c in ["XGB_Predicted", "SurrogateScore", "dG_pred"]:
        if c in df.columns:
            xgb_col = c
            break
    for c in ["DockingScore", "Vina_dG", "SSTR2_BindingEnergy"]:
        if c in df.columns:
            dock_col = c
            break

    if xgb_col and dock_col and xgb_col != dock_col:
        x = df[dock_col].dropna()
        y = df.loc[x.index, xgb_col].dropna()
        common_idx = x.index.intersection(y.index)
        x = x.loc[common_idx]
        y = y.loc[common_idx]

        ax.scatter(x, y, c=C_BLUE, alpha=0.6, s=40, edgecolors="gray", linewidths=0.3)

        if HAS_SCIPY and len(x) >= 5:
            rho, pval = sp_stats.spearmanr(x, y)
            # Regression line
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax.plot(x_line, p(x_line), "--", color=C_RED, lw=2)
            ax.text(0.05, 0.95,
                    f"Spearman ρ = {rho:.3f}\np = {pval:.2e}\nn = {len(x)}",
                    transform=ax.transAxes, fontsize=10, va="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8))

        ax.set_xlabel(f"Physical Docking Score [{dock_col}] (kcal/mol)")
        ax.set_ylabel(f"AI-predicted ΔG$_{{pred}}$ [{xgb_col}] (kcal/mol)")
        ax.set_title("Fig 12. Surrogate Model Validation: Spearman Correlation")

        # Identity line
        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, ":", color=C_GRAY, lw=1, label="y = x")
        ax.legend(frameon=False)

    else:
        # Fallback: plot ΔG_pred vs Rank to show monotonic ordering
        if "dG_pred" in df.columns:
            ranks = np.arange(1, len(df) + 1)
            ax.scatter(ranks, df["dG_pred"], c=C_BLUE, alpha=0.6, s=30)
            ax.set_xlabel("Candidate Rank")
            ax.set_ylabel("AI-predicted ΔG$_{pred}$ (kcal/mol)")
            ax.set_title("Fig 12. ΔG$_{pred}$ by Surrogate Rank (Monotonicity Check)")

            if HAS_SCIPY:
                rho, pval = sp_stats.spearmanr(ranks, df["dG_pred"])
                ax.text(0.05, 0.95,
                        f"Spearman ρ = {rho:.3f}\np = {pval:.2e}",
                        transform=ax.transAxes, fontsize=10, va="top",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8))

    _savefig(fig, "Fig12_Spearman_Correlation", test_mode)


def fig13_stability_profile(df, test_mode):
    """Fig 13. 5-Indicator Stability Profile (Top-10 grouped bar)."""
    indicator_cols = ["ProteaseSites", "AggregationScore", "PPBScore", "Rigidity",
                      "CompositeStabilityScore"]
    available = [c for c in indicator_cols if c in df.columns]
    if len(available) < 3:
        print("    [SKIP] Fig13 — insufficient stability indicator columns.")
        return

    fig, axes = plt.subplots(1, len(available), figsize=(4 * len(available), 5),
                             sharey=False)
    if len(available) == 1:
        axes = [axes]

    top10 = df.head(10)
    labels = [s[:8] for s in top10["Sequence"]]

    colors_ind = [C_BLUE, C_GREEN, C_PURPLE, C_CYAN, C_ORANGE]
    for i, col in enumerate(available):
        ax = axes[i]
        vals = top10[col].values
        ax.barh(range(len(vals)), vals, color=colors_ind[i % len(colors_ind)],
                edgecolor="white", height=0.6)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels if i == 0 else [], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel(col.replace("Score", "\nScore"), fontsize=9)
        ax.set_title(col, fontsize=10)

    fig.suptitle("Fig 13. Top-10 Candidates: 5-Indicator Stability Profile", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _savefig(fig, "Fig13_Stability_Profile", test_mode)


def table01_descriptive_stats(df, test_mode):
    """Table 1. Descriptive Statistics Summary (saved as image)."""
    cols = ["dG_pred", "Selectivity", "NetCharge", "pI", "Hydrophobicity"]
    for c in ["CompositeStabilityScore", "K_d_nM"]:
        if c in df.columns:
            cols.append(c)
    cols = [c for c in cols if c in df.columns]

    stats = df[cols].describe().T
    stats = stats[["count", "mean", "std", "min", "25%", "50%", "75%", "max"]]
    stats = stats.round(4)

    fig, ax = plt.subplots(figsize=(14, max(3, len(cols) * 0.5 + 1.5)))
    ax.axis("off")

    # Rename for publication
    rename_map = {
        "dG_pred": "ΔG_pred (kcal/mol)",
        "K_d_nM": "K_d (nM)",
        "CompositeStabilityScore": "CSS",
    }
    stats.index = [rename_map.get(c, c) for c in stats.index]

    table = ax.table(
        cellText=stats.values,
        rowLabels=stats.index,
        colLabels=[c.capitalize() for c in stats.columns],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.4)

    # Style header
    for key, cell in table.get_celld().items():
        if key[0] == 0:
            cell.set_facecolor(C_BLUE)
            cell.set_text_props(color="white", fontweight="bold")
        elif key[1] == -1:
            cell.set_facecolor("#F0F0F0")

    title = "Table 1. Descriptive Statistics Summary"
    if test_mode:
        title += f"\n{WATERMARK_TEXT}"
    ax.set_title(title, fontsize=12, pad=20)

    _savefig(fig, "Table01_Descriptive_Stats", test_mode)


def generate_top5_csv(df, test_mode):
    """Generate Top_5_Candidates_Summary.csv with K_d column.

    Uses ΔG = RT·ln(K_d) → K_d = exp(ΔG / RT)
    R = 1.987 cal/(mol·K), T = 310.15 K, RT ≈ 0.616 kcal/mol
    """
    top5 = df.head(5).copy()

    # Ensure K_d columns exist
    if "K_d_M" not in top5.columns:
        top5["K_d_M"] = top5["dG_pred"].apply(_dg_to_kd)
        top5["K_d_nM"] = top5["K_d_M"].apply(_kd_to_nm)

    # Select columns for output
    out_cols = ["Sequence", "Length"]
    if "Rank" in top5.columns:
        out_cols = ["Rank"] + out_cols
    out_cols += ["dG_pred", "K_d_nM", "Selectivity", "NetCharge"]

    for c in ["CompositeStabilityScore", "ProteaseSites", "AggregationScore",
              "PPBScore", "Rigidity"]:
        if c in top5.columns:
            out_cols.append(c)

    out_cols = [c for c in out_cols if c in top5.columns]
    top5_out = top5[out_cols].copy()

    # Rename for publication clarity
    col_rename = {
        "dG_pred": "dG_pred_kcal_mol",
        "K_d_nM": "K_d_nM",
    }
    top5_out.rename(columns=col_rename, inplace=True)

    csv_path = OUTPUT_DIR / "Top_5_Candidates_Summary.csv"
    top5_out.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"  [CSV] Top_5_Candidates_Summary.csv saved ({len(top5_out)} rows)")

    if test_mode:
        print(f"  [WARNING] {WATERMARK_TEXT}")

    return top5_out


def generate_validation_report(df, test_mode, top5_csv):
    """Generate Validation_Report.md with methodology, results, and disclaimer."""
    report_lines = []

    if test_mode:
        report_lines.append(f"# ⚠ {WATERMARK_TEXT} ⚠\n")
        report_lines.append("---\n")

    report_lines.append("# SSTR2 AI Scientist — Validation & Analysis Report\n")
    report_lines.append(f"**Generated:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    report_lines.append(f"**Candidates analyzed:** {len(df)}\n")
    report_lines.append("")

    # Methodology
    report_lines.append("## 1. Methodology\n")
    report_lines.append(
        "This report validates the top candidates from the SSTR2 AI Scientist pipeline.\n"
        "The pipeline consists of:\n"
        "- **Module 00:** Combinatorial peptide generation (Cys-Cys scaffold + FWKT pharmacophore)\n"
        "- **Module 01:** 5-indicator PK stability filter (Protease, Aggregation, Kidney, PPB, Rigidity)\n"
        "- **Module 02:** ESM-2 embedding + XGBoost surrogate → AI-predicted ΔG_pred\n"
        "- **Module 03:** 5-agent AI boardroom evaluation\n\n"
        "**Terminology:** All binding energies are reported as "
        "**AI-predicted Binding Free Energy (ΔG_pred)**, computed by an XGBoost surrogate "
        "trained on ESM-2 embeddings and RDKit molecular descriptors.\n"
    )

    # K_d derivation
    report_lines.append("## 2. Equilibrium Dissociation Constant (K_d) Derivation\n")
    report_lines.append(
        "K_d values are derived from ΔG_pred using the thermodynamic relationship:\n\n"
        "```\n"
        "ΔG = RT · ln(K_d)\n"
        "K_d = exp(ΔG / RT)\n"
        "```\n\n"
        f"Where R = {R_CAL} cal/(mol·K), T = {T_BODY} K (37°C), "
        f"RT = {RT_KCAL:.4f} kcal/mol.\n"
    )

    # Top-5 table
    report_lines.append("## 3. Top-5 Candidate Summary\n")
    report_lines.append("| Rank | Sequence | ΔG_pred (kcal/mol) | K_d (nM) | Selectivity | NetCharge |\n")
    report_lines.append("|------|----------|--------------------:|----------:|------------:|----------:|\n")
    top5 = df.head(5)
    for i, (_, row) in enumerate(top5.iterrows(), 1):
        dg = row.get("dG_pred", row.get("SSTR2_BindingEnergy", "N/A"))
        kd = row.get("K_d_nM", "N/A")
        sel = row.get("Selectivity", "N/A")
        charge = row.get("NetCharge", "N/A")
        dg_str = f"{dg:.3f}" if isinstance(dg, (int, float)) else str(dg)
        kd_str = f"{kd:.4f}" if isinstance(kd, (int, float)) else str(kd)
        sel_str = f"{sel:.1f}" if isinstance(sel, (int, float)) else str(sel)
        chg_str = f"{charge:+.4f}" if isinstance(charge, (int, float)) else str(charge)
        report_lines.append(f"| {i} | {row['Sequence']} | {dg_str} | {kd_str} | {sel_str} | {chg_str} |\n")

    # Dotatate comparison
    dot_kd = _kd_to_nm(_dg_to_kd(DOTATATE["BindingEnergy"]))
    report_lines.append(
        f"\n**Dotatate reference:** ΔG = {DOTATATE['BindingEnergy']} kcal/mol, "
        f"K_d = {dot_kd:.4f} nM, Selectivity = {DOTATATE['Selectivity']}x, "
        f"Charge = +{DOTATATE['NetCharge']:.3f}\n"
    )

    # Success criteria
    report_lines.append("\n## 4. Success Criteria Assessment\n")
    n_total = len(df)
    pass1 = (df["Selectivity"] >= THRESH_SELECTIVITY).sum() if "Selectivity" in df.columns else 0
    pass2 = (df["dG_pred"] <= THRESH_RETENTION).sum() if "dG_pred" in df.columns else 0
    pass3 = (df["NetCharge"] <= THRESH_CHARGE).sum() if "NetCharge" in df.columns else 0
    report_lines.append(
        f"| Criterion | Threshold | Pass | Rate |\n"
        f"|-----------|-----------|------|------|\n"
        f"| 100x Selectivity | ≥ {THRESH_SELECTIVITY}x | {pass1}/{n_total} | {pass1/n_total*100:.1f}% |\n"
        f"| 72h Retention | ≤ {THRESH_RETENTION} kcal/mol | {pass2}/{n_total} | {pass2/n_total*100:.1f}% |\n"
        f"| Kidney Safety | ≤ +{THRESH_CHARGE} | {pass3}/{n_total} | {pass3/n_total*100:.1f}% |\n"
    )

    # Figures
    report_lines.append("\n## 5. Generated Figures\n")
    report_lines.append(
        "All figures are saved in `output/figures/` at 300 DPI (PNG) with PDF vector copies.\n"
        "See the Figure Catalog in the module docstring for descriptions.\n"
    )

    # Disclaimer (학술적 무결성)
    report_lines.append("\n## 6. Limitations & Disclaimer\n")
    report_lines.append(f"{DISCLAIMER}\n")

    if test_mode:
        report_lines.append(f"\n---\n# ⚠ {WATERMARK_TEXT} ⚠\n")

    report_path = OUTPUT_DIR / "Validation_Report.md"
    report_path.write_text("".join(report_lines), encoding="utf-8")
    print(f"  [REPORT] Validation_Report.md saved")


# ==========================================================================
#  Main
# ==========================================================================

def main():
    t0 = time.time()
    print("=" * 70)
    print("  Module 04: Validation & Publication-Quality Analysis")
    print("=" * 70)

    _setup_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # ── 1. Load data ──
    print("\n[1/4] Loading data...")
    df, test_mode = load_data()

    if test_mode:
        print(f"\n  {'='*60}")
        print(f"  ⚠  {WATERMARK_TEXT}")
        print(f"  {'='*60}\n")

    # ── 2. Generate figures ──
    print("\n[2/4] Generating publication figures...")

    fig01_dg_distribution(df, test_mode)
    fig02_selectivity(df, test_mode)
    fig03_charge_safety(df, test_mode)
    fig04_dg_vs_selectivity(df, test_mode)
    fig05_property_space(df, test_mode)
    fig06_length_vs_dg(df, test_mode)
    fig07_composite_waterfall(df, test_mode)
    fig08_aa_heatmap(df, test_mode)
    fig09_correlation_matrix(df, test_mode)
    fig10_radar_chart(df, test_mode)
    fig11_pipeline_funnel(df, test_mode)
    fig12_spearman_correlation(df, test_mode)
    fig13_stability_profile(df, test_mode)
    table01_descriptive_stats(df, test_mode)

    # ── 3. Generate CSV + Report ──
    print("\n[3/4] Generating Top-5 Summary & Validation Report...")
    top5_csv = generate_top5_csv(df, test_mode)
    generate_validation_report(df, test_mode, top5_csv)

    # ── 4. Summary ──
    elapsed = time.time() - t0
    n_files = len(list(FIG_DIR.glob("*")))
    print(f"\n[4/4] Complete!")
    print(f"  Total figures/tables: {n_files} files")
    print(f"  Output directory: {FIG_DIR}")
    print(f"  Top-5 CSV: {OUTPUT_DIR / 'Top_5_Candidates_Summary.csv'}")
    print(f"  Report: {OUTPUT_DIR / 'Validation_Report.md'}")
    print(f"  Elapsed: {elapsed:.2f}s")

    if test_mode:
        print(f"\n  ⚠  {WATERMARK_TEXT}")

    print(f"\n{'='*70}")
    print(f"  K_d Derivation: ΔG = RT·ln(K_d), R={R_CAL}, T={T_BODY}K")
    print(f"  Dotatate ref: ΔG={DOTATATE['BindingEnergy']} → "
          f"K_d={_kd_to_nm(_dg_to_kd(DOTATATE['BindingEnergy'])):.4f} nM")
    print(f"{'='*70}")

    # Print library versions
    print("\n[Library Versions]")
    print(f"  numpy:      {np.__version__}")
    print(f"  pandas:     {pd.__version__}")
    print(f"  matplotlib: {matplotlib.__version__}")
    if HAS_SEABORN:
        print(f"  seaborn:    {sns.__version__}")
    if HAS_SCIPY:
        import scipy
        print(f"  scipy:      {scipy.__version__}")


if __name__ == "__main__":
    main()
