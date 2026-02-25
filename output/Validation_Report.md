# SSTR2 AI Scientist — Validation & Analysis Report
**Generated:** 2026-02-25 17:34:53
**Candidates analyzed:** 100
## 1. Methodology
This report validates the top candidates from the SSTR2 AI Scientist pipeline.
The pipeline consists of:
- **Module 00:** Combinatorial peptide generation (Cys-Cys scaffold + FWKT pharmacophore)
- **Module 01:** 5-indicator PK stability filter (Protease, Aggregation, Kidney, PPB, Rigidity)
- **Module 02:** ESM-2 embedding + XGBoost surrogate → AI-predicted ΔG_pred
- **Module 03:** 5-agent AI boardroom evaluation

**Terminology:** All binding energies are reported as **AI-predicted Binding Free Energy (ΔG_pred)**, computed by an XGBoost surrogate trained on ESM-2 embeddings and RDKit molecular descriptors.
## 2. Equilibrium Dissociation Constant (K_d) Derivation
K_d values are derived from ΔG_pred using the thermodynamic relationship:

```
ΔG = RT · ln(K_d)
K_d = exp(ΔG / RT)
```

Where R = 1.987 cal/(mol·K), T = 310.15 K (37°C), RT = 0.6163 kcal/mol.
## 3. Top-5 Candidate Summary
| Rank | Sequence | ΔG_pred (kcal/mol) | K_d (nM) | Selectivity | NetCharge |
|------|----------|--------------------:|----------:|------------:|----------:|
| 1 | CAFWKTCE | -13.410 | 0.3546 | 30828.1 | -0.2898 |
| 2 | CFWKTPEC | -13.330 | 0.4038 | 28975.6 | -0.2898 |
| 3 | CPFWKTEC | -13.058 | 0.6278 | 23501.4 | -0.2898 |
| 4 | CFWKTDCW | -12.969 | 0.7253 | 21948.6 | -0.2903 |
| 5 | PPCDFWKTCF | -12.605 | 1.3093 | 16578.8 | -0.2903 |

**Dotatate reference:** ΔG = -8.5 kcal/mol, K_d = 1023.0733 nM, Selectivity = 50.0x, Charge = +0.707

## 4. Success Criteria Assessment
| Criterion | Threshold | Pass | Rate |
|-----------|-----------|------|------|
| 100x Selectivity | ≥ 100.0x | 100/100 | 100.0% |
| 72h Retention | ≤ -7.0 kcal/mol | 100/100 | 100.0% |
| Kidney Safety | ≤ +0.707 | 100/100 | 100.0% |

## 5. Generated Figures
All figures are saved in `output/figures/` at 300 DPI (PNG) with PDF vector copies.
See the Figure Catalog in the module docstring for descriptions.

## 6. Limitations & Disclaimer
**Disclaimer:** This Surrogate Model combines ESM-2 embeddings, RDKit molecular descriptors, and XGBoost regression. The reported ΔG_pred values are AI-predicted approximations, NOT experimentally measured binding affinities. All final candidates MUST be validated through physical molecular docking (e.g., AutoDock Vina, Glide) and molecular dynamics (MD) simulations before any experimental synthesis or in vivo testing.
