#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
파이프라인 시뮬레이션: 현실적인 출력 데이터 생성 + Module 04 실행.
RDKit 없이도 동작하도록 peptides 라이브러리만 사용.
"""

import sys
import math
import time
import numpy as np
import pandas as pd
from pathlib import Path

if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)
(OUTPUT_DIR / "figures").mkdir(exist_ok=True)

rng = np.random.RandomState(42)

VARIABLE_AA = list("ADEFGHILMNPQSTVWY")
RT_KCAL = 1.987 * 310.15 / 1000.0


def generate_sequence(rng, length):
    n_var = length - 6
    dividers = sorted(rng.randint(0, n_var + 1, size=3))
    parts = [dividers[0], dividers[1] - dividers[0],
             dividers[2] - dividers[1], n_var - dividers[2]]
    n_tail = "".join(rng.choice(VARIABLE_AA, size=parts[0]))
    inner_l = "".join(rng.choice(VARIABLE_AA, size=parts[1]))
    inner_r = "".join(rng.choice(VARIABLE_AA, size=parts[2]))
    c_tail = "".join(rng.choice(VARIABLE_AA, size=parts[3]))
    return n_tail + "C" + inner_l + "FWKT" + inner_r + "C" + c_tail


# ==========================================================================
print("=" * 60)
print("  Pipeline Simulation (No RDKit Required)")
print("=" * 60)

# ── Module 00 ─────────────────────────────────────────────────────────────
print("\n>>> Module 00: Peptide Library Generation")
t0 = time.time()

sequences = set()
while len(sequences) < 10000:
    length = rng.randint(8, 15)
    seq = generate_sequence(rng, length)
    sequences.add(seq)

sequences = list(sequences)
rng.shuffle(sequences)

df_raw = pd.DataFrame({"Sequence": sequences, "Length": [len(s) for s in sequences]})
df_raw.to_csv(OUTPUT_DIR / "SSTR2_Raw_Library.csv", index=False)
print(f"  Generated: {len(df_raw):,} sequences ({time.time()-t0:.1f}s)")

# ── Module 01 ─────────────────────────────────────────────────────────────
print("\n>>> Module 01: PK Stability Filtering")

try:
    from peptides import Peptide
    HAS_PEPTIDES = True
    print("  peptides library: OK")
except ImportError:
    HAS_PEPTIDES = False
    print("  [WARNING] peptides not installed, using approximation")

results = []
for idx, seq in enumerate(sequences):
    if idx % 2000 == 0:
        print(f"\r  Processing: {idx:,}/{len(sequences):,}", end="", flush=True)

    # Protease sites
    protease_sites = 0
    fwkt_start = seq.find("FWKT")
    for i in range(len(seq) - 1):
        aa = seq[i]
        if aa in "KR" and i != (fwkt_start + 2 if fwkt_start >= 0 else -1):
            protease_sites += 1
        if aa in "FWYL" and not (fwkt_start >= 0 and fwkt_start <= i <= fwkt_start + 3):
            protease_sites += 1

    if HAS_PEPTIDES:
        pep = Peptide(seq)
        charge = pep.charge(pH=7.4)
        pi_val = pep.isoelectric_point()
        hydro = pep.hydrophobicity("Eisenberg")
        boman = pep.boman()
        kd_hydro = pep.hydrophobicity("KyteDoolittle")
    else:
        # Approximation without peptides library
        pos_aa = sum(1 for a in seq if a in "KRH")
        neg_aa = sum(1 for a in seq if a in "DE")
        charge = round(pos_aa * 0.8 - neg_aa * 0.9 + rng.normal(0, 0.1), 4)
        pi_val = round(6.0 + pos_aa * 0.5 - neg_aa * 0.4 + rng.normal(0, 0.3), 4)
        hydro = round(sum(-0.3 if a in "DEKRH" else 0.2 if a in "FWYILV" else 0.0
                         for a in seq) / len(seq) + rng.normal(0, 0.05), 4)
        boman = round(1.0 + rng.normal(0, 0.5), 4)
        kd_hydro = round(hydro * 1.2 + rng.normal(0, 0.1), 4)

    agg_score = round(max(boman, 0) + max(kd_hydro, 0), 4)

    # PPB (simplified without RDKit)
    n_polar = sum(1 for a in seq if a in "DEKRHNSTYQC")
    logp_approx = kd_hydro * len(seq) if HAS_PEPTIDES else hydro * len(seq)
    tpsa_approx = n_polar * 40.0 + 60.0
    ppb = round(logp_approx - 0.01 * tpsa_approx, 4)

    # Rigidity (simplified)
    rigidity = max(0, (len(seq) - 1) * 2 - seq.count("P") - (seq.count("C") // 2) * 4)

    results.append({
        "Sequence": seq, "Length": len(seq),
        "NetCharge": round(charge, 4), "pI": round(pi_val, 4),
        "Hydrophobicity": round(hydro, 4),
        "ProteaseSites": protease_sites,
        "AggregationScore": agg_score,
        "PPBScore": ppb,
        "Rigidity": int(rigidity),
    })

print(f"\r  Processed: {len(sequences):,} sequences")

df_all = pd.DataFrame(results)

# Hard cutoffs
mask = (df_all["NetCharge"] <= 1.0) & (df_all["pI"] <= 9.3)
df_pass = df_all[mask].copy()
print(f"  Hard cutoff pass: {len(df_pass):,} / {len(df_all):,}")

# CSS
for col in ["ProteaseSites", "AggregationScore", "PPBScore", "Rigidity"]:
    mn, mx = df_pass[col].min(), df_pass[col].max()
    s_name = "S_" + col.replace("Score", "").replace("Sites", "")
    if mx > mn:
        df_pass[s_name] = (1.0 - (df_pass[col] - mn) / (mx - mn)).round(4)
    else:
        df_pass[s_name] = 1.0

kidney_penalty = np.maximum(df_pass["NetCharge"].values, 0.0)
mn, mx = kidney_penalty.min(), kidney_penalty.max()
if mx > mn:
    df_pass["S_Kidney"] = (1.0 - (kidney_penalty - mn) / (mx - mn)).round(4)
else:
    df_pass["S_Kidney"] = 1.0

df_pass["CompositeStabilityScore"] = (
    0.20 * df_pass["S_Protease"] +
    0.15 * df_pass["S_Aggregation"] +
    0.30 * df_pass["S_Kidney"] +
    0.20 * df_pass["S_PPB"] +
    0.15 * df_pass["S_Rigidity"]
).round(4)

df_top100 = df_pass.nlargest(100, "CompositeStabilityScore").copy()
df_top100.insert(0, "Rank", range(1, 101))
df_top100.to_csv(OUTPUT_DIR / "SSTR2_Filtered_Top100.csv", index=False)
print(f"  Top-100 CSS: [{df_top100['CompositeStabilityScore'].min():.4f}, "
      f"{df_top100['CompositeStabilityScore'].max():.4f}]")
print(f"  Mean charge: {df_top100['NetCharge'].mean():+.4f}")

# ── Module 02 ─────────────────────────────────────────────────────────────
print("\n>>> Module 02: Docking + Surrogate Model Simulation")

seqs_100 = df_top100["Sequence"].tolist()

def simulate_docking(seq, rng):
    score = -6.5
    c_pos = [i for i, aa in enumerate(seq) if aa == "C"]
    if len(c_pos) >= 2:
        ring_size = c_pos[1] - c_pos[0] - 1
        if 5 <= ring_size <= 7:
            score -= 2.5
        elif 4 <= ring_size <= 8:
            score -= 1.5
        ring_seq = seq[c_pos[0]+1:c_pos[1]]
        score -= sum(1 for a in ring_seq if a in "FWYILV") * 0.25
        score += sum(1 for a in ring_seq if a in "DE") * 0.3
        if "P" in ring_seq:
            score -= 0.5
    fwkt_idx = seq.find("FWKT")
    if fwkt_idx >= 0:
        center_dist = abs((fwkt_idx + 2) / len(seq) - 0.5)
        score -= (0.5 - center_dist) * 0.8
    if len(seq) > 12:
        score += 0.2 * (len(seq) - 12)
    if seq[0] in "AGSTP":
        score -= 0.3
    score += rng.normal(0, 0.35)
    return round(max(score, -12.5), 3)

docking_scores = np.array([simulate_docking(s, rng) for s in seqs_100])

def simulate_selectivity(dg, charge):
    factors = [0.55, 0.50, 0.35, 0.70]
    nt_scores = [min(dg * f, -0.5) for f in factors]
    ddg = np.mean(nt_scores) - dg
    sel = np.exp(ddg / RT_KCAL)
    if charge < 0:
        sel *= 1.1
    return round(sel, 1)

selectivities = np.array([
    simulate_selectivity(docking_scores[i], df_top100.iloc[i]["NetCharge"])
    for i in range(100)
])

log_sel = np.log10(np.clip(selectivities, 1.0, None))
composite = np.round(docking_scores - log_sel, 3)

df_m02 = df_top100.copy()
df_m02["SSTR2_BindingEnergy"] = docking_scores
df_m02["Selectivity"] = selectivities
df_m02["CompositeScore"] = composite
df_m02 = df_m02.sort_values("CompositeScore").reset_index(drop=True)
df_m02["Rank"] = range(1, 101)

col_order = ["Rank", "Sequence", "Length", "SSTR2_BindingEnergy", "Selectivity",
             "CompositeScore", "NetCharge", "pI", "Hydrophobicity",
             "ProteaseSites", "AggregationScore", "PPBScore", "Rigidity",
             "CompositeStabilityScore"]
col_order = [c for c in col_order if c in df_m02.columns]
df_m02 = df_m02[col_order]
df_m02.to_csv(OUTPUT_DIR / "SSTR2_Top100_Boardroom_Ready.csv", index=False)

print(f"  dG range: [{docking_scores.min():.3f}, {docking_scores.max():.3f}]")
print(f"  Selectivity range: [{selectivities.min():.1f}x, {selectivities.max():.1f}x]")
print(f"  Pass dG <= -7.0: {(docking_scores <= -7.0).sum()}/100")
print(f"  Pass Sel >= 100x: {(selectivities >= 100).sum()}/100")

print(f"\n  [Top 5]")
for _, row in df_m02.head(5).iterrows():
    print(f"    #{int(row['Rank']):>2} {row['Sequence']:<16} dG={row['SSTR2_BindingEnergy']:>7.3f}  "
          f"Sel={row['Selectivity']:>7.1f}x  Chg={row['NetCharge']:>+7.4f}")

# ── Module 03 ─────────────────────────────────────────────────────────────
print("\n>>> Module 03: AI Boardroom Simulation")

top5 = df_m02.head(5)
approved = []
transcript_lines = []
transcript_lines.append("=" * 70)
transcript_lines.append("  SSTR2 AI Scientist - Boardroom Transcript")
transcript_lines.append(f"  Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
transcript_lines.append("=" * 70 + "\n")

for i, (_, cand) in enumerate(top5.iterrows()):
    seq = cand["Sequence"]
    dg = cand["SSTR2_BindingEnergy"]
    sel = cand["Selectivity"]
    charge = cand["NetCharge"]
    css = cand["CompositeStabilityScore"]
    c_pos = [j for j, aa in enumerate(seq) if aa == "C"]
    ring_size = c_pos[1] - c_pos[0] - 1 if len(c_pos) >= 2 else 0

    transcript_lines.append("-" * 50)
    transcript_lines.append(f"[Evidence Curator] (Seq: {seq})")
    transcript_lines.append(f"  Structure: {len(seq)}AA, Cys{c_pos[0]+1}-Cys{c_pos[1]+1}, ring={ring_size}")
    transcript_lines.append(f"  FWKT at pos {seq.find('FWKT')+1}, PDB 7T10 compatible")
    transcript_lines.append(f"  Trp→TM3-TM7 pocket, Lys→Asp3.32 salt bridge")

    transcript_lines.append("-" * 50)
    struct_pass = dg <= -7.0 and sel >= 100.0
    verdict_s = "PASS" if struct_pass else "REJECT"
    transcript_lines.append(f"[Structure Analyst] (Seq: {seq}) [{verdict_s}]")
    transcript_lines.append(f"  Selectivity: {sel:.1f}x {'>=100' if sel>=100 else '<100'}")
    transcript_lines.append(f"  Retention: {dg:.3f} {'<=-7.0' if dg<=-7.0 else '>-7.0'}")
    transcript_lines.append(f"  VERDICT: {verdict_s}")

    transcript_lines.append("-" * 50)
    pk_pass = charge <= 0.8 and css >= 0.3
    verdict_p = "PASS" if pk_pass else "REJECT"
    transcript_lines.append(f"[PK Expert] (Seq: {seq}) [{verdict_p}]")
    transcript_lines.append(f"  Charge: {charge:+.4f}, CSS: {css:.4f}")
    transcript_lines.append(f"  VERDICT: {verdict_p}")

    if struct_pass and pk_pass:
        approved.append(cand.to_dict())
        transcript_lines.append("-" * 50)
        transcript_lines.append(f"[Experiment Planner] (Seq: {seq}) [APPROVED]")
        transcript_lines.append(f"  Protocol: DOTA-NHS → 177Lu label → AR42J binding → xenograft biodist")
    else:
        transcript_lines.append("-" * 50)
        transcript_lines.append(f"[System] (Seq: {seq}) [REJECT]")
        transcript_lines.append(f"  Reason: {'low selectivity/retention' if not struct_pass else 'PK concerns'}")

    transcript_lines.append("")

transcript_lines.append("=" * 70)
transcript_lines.append(f"  Summary: {len(approved)}/{len(top5)} approved")
transcript_lines.append("=" * 70)

(OUTPUT_DIR / "Boardroom_Transcript.log").write_text("\n".join(transcript_lines), encoding="utf-8")
print(f"  Approved: {len(approved)}/{len(top5)}")

# Final Report
report = []
report.append("# SSTR2 AI Scientist — Final Winning Peptides Report\n\n")
report.append(f"**Generated:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
report.append("## Approved Candidates\n\n")
report.append("| Sequence | dG_pred | Selectivity | NetCharge | CSS |\n")
report.append("|----------|--------:|------------:|----------:|----:|\n")
for a in approved:
    report.append(f"| {a['Sequence']} | {a['SSTR2_BindingEnergy']:.3f} | "
                  f"{a['Selectivity']:.1f}x | {a['NetCharge']:+.4f} | "
                  f"{a['CompositeStabilityScore']:.4f} |\n")
report.append("\n## Selection Rationale\n\n")
for i, a in enumerate(approved, 1):
    seq = a["Sequence"]
    kd_nm = math.exp(a["SSTR2_BindingEnergy"] / RT_KCAL) * 1e9
    report.append(f"### Candidate {i}: {seq}\n")
    report.append(f"- dG = {a['SSTR2_BindingEnergy']:.3f} kcal/mol, Kd = {kd_nm:.2f} nM\n")
    report.append(f"- Selectivity = {a['Selectivity']:.1f}x, Charge = {a['NetCharge']:+.4f}\n\n")

(OUTPUT_DIR / "Final_Winning_Peptides_Report.md").write_text("".join(report), encoding="utf-8")
print(f"  Saved: Final_Winning_Peptides_Report.md")

# ── Module 04 ─────────────────────────────────────────────────────────────
print("\n>>> Module 04: Validation & Figure Generation")
print("  Running 04_validate_and_analyze.py ...")

import subprocess
result = subprocess.run(
    [sys.executable, "04_validate_and_analyze.py"],
    capture_output=True, text=True, timeout=300
)
print(result.stdout)
if result.stderr:
    print("  [STDERR]:", result.stderr[-500:])

if result.returncode != 0:
    print(f"  [WARNING] Module 04 exit code: {result.returncode}")
    print("  Figures may be incomplete, but CSV data is ready for dashboard.")

# ── Summary ───────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("  SIMULATION COMPLETE")
print("=" * 60)
for f in sorted(OUTPUT_DIR.glob("*")):
    if f.is_file():
        print(f"  {f.name:45s} {f.stat().st_size/1024:>7.1f} KB")
fig_count = len(list((OUTPUT_DIR / "figures").glob("*")))
print(f"\n  Figures: {fig_count} files")
print(f"\n  Next: streamlit run dashboard.py")
print("=" * 60)
