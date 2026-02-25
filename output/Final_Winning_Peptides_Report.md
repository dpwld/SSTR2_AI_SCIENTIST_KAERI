# SSTR2 AI Scientist — Final Winning Peptides Report

**Generated:** 2026-02-25 17:34:27

## 1. Project Overview

This report presents the final candidates selected by the 5-agent AI Boardroom system for next-generation SSTR2-targeted Peptide Receptor Radionuclide Therapy (PRRT). Each candidate has been evaluated against three mandatory success criteria by specialized AI agents (Evidence Curator, Structure Analyst, PK Expert, Design Generator, Experiment Planner).

## 2. Success Criteria

| # | Criterion | Threshold | Rationale |
|---|-----------|-----------|----------|
| 1 | **Target Selectivity** | Selectivity >= 100x vs SSTR1/3/4/5 | Minimize off-target radiation to healthy tissue |
| 2 | **72h Retention** | Binding energy <= -7.0 kcal/mol | k_off < 3.86e-6 /s for sustained Lu-177 irradiation |
| 3 | **Kidney Safety** | Net charge <= Dotatate (+0.707) | Reduce megalin-mediated renal tubular reabsorption |

## 3. Approved Candidates (5)

| # | Sequence | Binding Energy (kcal/mol) | Selectivity | Net Charge | Composite Score |
|---|----------|--------------------------|-------------|------------|----------------|
| 1 | `CAFWKTCE` | -13.41 | 30828.1x | -0.2898 | -17.899 |
| 2 | `CFWKTPEC` | -13.33 | 28975.6x | -0.2898 | -17.792 |
| 3 | `CPFWKTEC` | -13.058 | 23501.4x | -0.2898 | -17.429 |
| 4 | `CFWKTDCW` | -12.969 | 21948.6x | -0.2903 | -17.311 |
| 5 | `PPCDFWKTCF` | -12.605 | 16578.8x | -0.2903 | -16.825 |

## 4. Selection Rationale

Each candidate was approved based on the following structural, biophysical, and pharmacokinetic evidence:

### Candidate #1: `CAFWKTCE`

**Structural Basis:**
- Scaffold: 8 AA, disulfide ring Cys1-Cys7 (ring size: 5 residues)
- FWKT pharmacophore at position 3-6: Trp inserts into SSTR2 TM3-TM7 hydrophobic pocket (PDB 7T10); Lys forms salt bridge with Asp3.32
- Beta-turn constraint from disulfide ring reduces conformational entropy penalty upon binding

**Criterion 1 — 100x Selectivity:** Selectivity = 30828.1x (threshold: >= 100x). The FWKT motif provides SSTR2-specific interactions that are poorly accommodated by SSTR1 (narrow Trp pocket), SSTR3 (ECL2 steric clash), SSTR4 (divergent binding site), and SSTR5 (weakened Lys salt bridge).

**Criterion 2 — 72h Retention:** Binding energy = -13.41 kcal/mol (threshold: <= -7.0). Strong binding yields k_off < 3.86e-6 /s, enabling >72h target residence for effective Lu-177 radiation delivery.

**Criterion 3 — Kidney Safety:** Net charge (-0.2898) is negative — substantially improved over Dotatate (+0.707). Reduced megalin-mediated renal reabsorption expected.

**Composite Score:** -17.899 (binding energy - log10(selectivity); lower = better)

---

### Candidate #2: `CFWKTPEC`

**Structural Basis:**
- Scaffold: 8 AA, disulfide ring Cys1-Cys8 (ring size: 6 residues)
- FWKT pharmacophore at position 2-5: Trp inserts into SSTR2 TM3-TM7 hydrophobic pocket (PDB 7T10); Lys forms salt bridge with Asp3.32
- Beta-turn constraint from disulfide ring reduces conformational entropy penalty upon binding

**Criterion 1 — 100x Selectivity:** Selectivity = 28975.6x (threshold: >= 100x). The FWKT motif provides SSTR2-specific interactions that are poorly accommodated by SSTR1 (narrow Trp pocket), SSTR3 (ECL2 steric clash), SSTR4 (divergent binding site), and SSTR5 (weakened Lys salt bridge).

**Criterion 2 — 72h Retention:** Binding energy = -13.33 kcal/mol (threshold: <= -7.0). Strong binding yields k_off < 3.86e-6 /s, enabling >72h target residence for effective Lu-177 radiation delivery.

**Criterion 3 — Kidney Safety:** Net charge (-0.2898) is negative — substantially improved over Dotatate (+0.707). Reduced megalin-mediated renal reabsorption expected.

**Composite Score:** -17.792 (binding energy - log10(selectivity); lower = better)

---

### Candidate #3: `CPFWKTEC`

**Structural Basis:**
- Scaffold: 8 AA, disulfide ring Cys1-Cys8 (ring size: 6 residues)
- FWKT pharmacophore at position 3-6: Trp inserts into SSTR2 TM3-TM7 hydrophobic pocket (PDB 7T10); Lys forms salt bridge with Asp3.32
- Beta-turn constraint from disulfide ring reduces conformational entropy penalty upon binding

**Criterion 1 — 100x Selectivity:** Selectivity = 23501.4x (threshold: >= 100x). The FWKT motif provides SSTR2-specific interactions that are poorly accommodated by SSTR1 (narrow Trp pocket), SSTR3 (ECL2 steric clash), SSTR4 (divergent binding site), and SSTR5 (weakened Lys salt bridge).

**Criterion 2 — 72h Retention:** Binding energy = -13.058 kcal/mol (threshold: <= -7.0). Strong binding yields k_off < 3.86e-6 /s, enabling >72h target residence for effective Lu-177 radiation delivery.

**Criterion 3 — Kidney Safety:** Net charge (-0.2898) is negative — substantially improved over Dotatate (+0.707). Reduced megalin-mediated renal reabsorption expected.

**Composite Score:** -17.429 (binding energy - log10(selectivity); lower = better)

---

### Candidate #4: `CFWKTDCW`

**Structural Basis:**
- Scaffold: 8 AA, disulfide ring Cys1-Cys7 (ring size: 5 residues)
- FWKT pharmacophore at position 2-5: Trp inserts into SSTR2 TM3-TM7 hydrophobic pocket (PDB 7T10); Lys forms salt bridge with Asp3.32
- Beta-turn constraint from disulfide ring reduces conformational entropy penalty upon binding

**Criterion 1 — 100x Selectivity:** Selectivity = 21948.6x (threshold: >= 100x). The FWKT motif provides SSTR2-specific interactions that are poorly accommodated by SSTR1 (narrow Trp pocket), SSTR3 (ECL2 steric clash), SSTR4 (divergent binding site), and SSTR5 (weakened Lys salt bridge).

**Criterion 2 — 72h Retention:** Binding energy = -12.969 kcal/mol (threshold: <= -7.0). Strong binding yields k_off < 3.86e-6 /s, enabling >72h target residence for effective Lu-177 radiation delivery.

**Criterion 3 — Kidney Safety:** Net charge (-0.2903) is negative — substantially improved over Dotatate (+0.707). Reduced megalin-mediated renal reabsorption expected.

**Composite Score:** -17.311 (binding energy - log10(selectivity); lower = better)

---

### Candidate #5: `PPCDFWKTCF`

**Structural Basis:**
- Scaffold: 10 AA, disulfide ring Cys3-Cys9 (ring size: 5 residues)
- FWKT pharmacophore at position 5-8: Trp inserts into SSTR2 TM3-TM7 hydrophobic pocket (PDB 7T10); Lys forms salt bridge with Asp3.32
- Beta-turn constraint from disulfide ring reduces conformational entropy penalty upon binding

**Criterion 1 — 100x Selectivity:** Selectivity = 16578.8x (threshold: >= 100x). The FWKT motif provides SSTR2-specific interactions that are poorly accommodated by SSTR1 (narrow Trp pocket), SSTR3 (ECL2 steric clash), SSTR4 (divergent binding site), and SSTR5 (weakened Lys salt bridge).

**Criterion 2 — 72h Retention:** Binding energy = -12.605 kcal/mol (threshold: <= -7.0). Strong binding yields k_off < 3.86e-6 /s, enabling >72h target residence for effective Lu-177 radiation delivery.

**Criterion 3 — Kidney Safety:** Net charge (-0.2903) is negative — substantially improved over Dotatate (+0.707). Reduced megalin-mediated renal reabsorption expected.

**Composite Score:** -16.825 (binding energy - log10(selectivity); lower = better)

---

## 5. Preclinical Experiment Plans

### Candidate #1: `CAFWKTCE`

```
[Experiment Planner - Preclinical Protocol]
Candidate: CAFWKTCE | Binding: -13.41 kcal/mol | Sel: 30828.1x | Charge: -0.2898

1. DOTA-NHS Coupling:
   - Conjugate DOTA-NHS ester to N-terminal alpha-amine
   - Conditions: pH 8.5 carbonate buffer, 25C, 2h
   - Purify by RP-HPLC (C18), confirm purity >95% by LC-MS

2. Lu-177 Radiolabeling:
   - DOTA-peptide + [Lu-177]LuCl3, 0.4M acetate buffer pH 5.0, 95C, 30min
   - Radiochemical purity >98% by iTLC-SG

3. In Vitro Binding:
   - AR42J cells (endogenous SSTR2): saturation binding for Kd
   - Selectivity panel: competition binding on CHO-SSTR1/3/4/5

4. In Vivo Biodistribution (AR42J xenograft nude mice, n=5/timepoint):
   - Timepoints: 1h, 4h, 24h, 48h, 72h post-injection
   - Measure %ID/g: tumor, kidneys, liver, blood, spleen
   - Calculate tumor-to-kidney ratio at each timepoint

5. Success Criteria:
   - Tumor uptake: >10 %ID/g at 4h
   - 72h retention: tumor %ID/g decrease <50% vs 4h
   - Tumor/kidney ratio: >1.2x vs [Lu-177]DOTATATE (20% improvement)
```

### Candidate #2: `CFWKTPEC`

```
[Experiment Planner - Preclinical Protocol]
Candidate: CFWKTPEC | Binding: -13.33 kcal/mol | Sel: 28975.6x | Charge: -0.2898

1. DOTA-NHS Coupling:
   - Conjugate DOTA-NHS ester to N-terminal alpha-amine
   - Conditions: pH 8.5 carbonate buffer, 25C, 2h
   - Purify by RP-HPLC (C18), confirm purity >95% by LC-MS

2. Lu-177 Radiolabeling:
   - DOTA-peptide + [Lu-177]LuCl3, 0.4M acetate buffer pH 5.0, 95C, 30min
   - Radiochemical purity >98% by iTLC-SG

3. In Vitro Binding:
   - AR42J cells (endogenous SSTR2): saturation binding for Kd
   - Selectivity panel: competition binding on CHO-SSTR1/3/4/5

4. In Vivo Biodistribution (AR42J xenograft nude mice, n=5/timepoint):
   - Timepoints: 1h, 4h, 24h, 48h, 72h post-injection
   - Measure %ID/g: tumor, kidneys, liver, blood, spleen
   - Calculate tumor-to-kidney ratio at each timepoint

5. Success Criteria:
   - Tumor uptake: >10 %ID/g at 4h
   - 72h retention: tumor %ID/g decrease <50% vs 4h
   - Tumor/kidney ratio: >1.2x vs [Lu-177]DOTATATE (20% improvement)
```

### Candidate #3: `CPFWKTEC`

```
[Experiment Planner - Preclinical Protocol]
Candidate: CPFWKTEC | Binding: -13.058 kcal/mol | Sel: 23501.4x | Charge: -0.2898

1. DOTA-NHS Coupling:
   - Conjugate DOTA-NHS ester to N-terminal alpha-amine
   - Conditions: pH 8.5 carbonate buffer, 25C, 2h
   - Purify by RP-HPLC (C18), confirm purity >95% by LC-MS

2. Lu-177 Radiolabeling:
   - DOTA-peptide + [Lu-177]LuCl3, 0.4M acetate buffer pH 5.0, 95C, 30min
   - Radiochemical purity >98% by iTLC-SG

3. In Vitro Binding:
   - AR42J cells (endogenous SSTR2): saturation binding for Kd
   - Selectivity panel: competition binding on CHO-SSTR1/3/4/5

4. In Vivo Biodistribution (AR42J xenograft nude mice, n=5/timepoint):
   - Timepoints: 1h, 4h, 24h, 48h, 72h post-injection
   - Measure %ID/g: tumor, kidneys, liver, blood, spleen
   - Calculate tumor-to-kidney ratio at each timepoint

5. Success Criteria:
   - Tumor uptake: >10 %ID/g at 4h
   - 72h retention: tumor %ID/g decrease <50% vs 4h
   - Tumor/kidney ratio: >1.2x vs [Lu-177]DOTATATE (20% improvement)
```

### Candidate #4: `CFWKTDCW`

```
[Experiment Planner - Preclinical Protocol]
Candidate: CFWKTDCW | Binding: -12.969 kcal/mol | Sel: 21948.6x | Charge: -0.2903

1. DOTA-NHS Coupling:
   - Conjugate DOTA-NHS ester to N-terminal alpha-amine
   - Conditions: pH 8.5 carbonate buffer, 25C, 2h
   - Purify by RP-HPLC (C18), confirm purity >95% by LC-MS

2. Lu-177 Radiolabeling:
   - DOTA-peptide + [Lu-177]LuCl3, 0.4M acetate buffer pH 5.0, 95C, 30min
   - Radiochemical purity >98% by iTLC-SG

3. In Vitro Binding:
   - AR42J cells (endogenous SSTR2): saturation binding for Kd
   - Selectivity panel: competition binding on CHO-SSTR1/3/4/5

4. In Vivo Biodistribution (AR42J xenograft nude mice, n=5/timepoint):
   - Timepoints: 1h, 4h, 24h, 48h, 72h post-injection
   - Measure %ID/g: tumor, kidneys, liver, blood, spleen
   - Calculate tumor-to-kidney ratio at each timepoint

5. Success Criteria:
   - Tumor uptake: >10 %ID/g at 4h
   - 72h retention: tumor %ID/g decrease <50% vs 4h
   - Tumor/kidney ratio: >1.2x vs [Lu-177]DOTATATE (20% improvement)
```

### Candidate #5: `PPCDFWKTCF`

```
[Experiment Planner - Preclinical Protocol]
Candidate: PPCDFWKTCF | Binding: -12.605 kcal/mol | Sel: 16578.8x | Charge: -0.2903

1. DOTA-NHS Coupling:
   - Conjugate DOTA-NHS ester to N-terminal alpha-amine
   - Conditions: pH 8.5 carbonate buffer, 25C, 2h
   - Purify by RP-HPLC (C18), confirm purity >95% by LC-MS

2. Lu-177 Radiolabeling:
   - DOTA-peptide + [Lu-177]LuCl3, 0.4M acetate buffer pH 5.0, 95C, 30min
   - Radiochemical purity >98% by iTLC-SG

3. In Vitro Binding:
   - AR42J cells (endogenous SSTR2): saturation binding for Kd
   - Selectivity panel: competition binding on CHO-SSTR1/3/4/5

4. In Vivo Biodistribution (AR42J xenograft nude mice, n=5/timepoint):
   - Timepoints: 1h, 4h, 24h, 48h, 72h post-injection
   - Measure %ID/g: tumor, kidneys, liver, blood, spleen
   - Calculate tumor-to-kidney ratio at each timepoint

5. Success Criteria:
   - Tumor uptake: >10 %ID/g at 4h
   - 72h retention: tumor %ID/g decrease <50% vs 4h
   - Tumor/kidney ratio: >1.2x vs [Lu-177]DOTATATE (20% improvement)
```

## 6. Agent Evaluation Summary

| Agent | PASS | REJECT |
|-------|------|--------|
| Structure Analyst | 5 | 0 |
| PK Expert | 5 | 0 |

---
*Generated by SSTR2 AI Scientist Boardroom System*
