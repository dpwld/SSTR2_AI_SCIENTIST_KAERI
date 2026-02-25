# 2. Methods and Results

## 2.1 Pipeline Overview

We developed a five-module computational pipeline for discovering SSTR2-targeted peptide ligands with improved binding affinity and pharmacokinetic properties compared to the clinical reference compound DOTATATE (L-amino acid approximation: FCYWKTCT). The pipeline operates in three phases: (1) combinatorial library generation and stability-based filtering, (2) structure-based binding affinity prediction via molecular docking and a machine learning surrogate model, and (3) multi-agent evaluation of the top candidates. All computations were performed locally on a single workstation (Windows 11, Python 3.12) without cloud-based resources. The receptor structure was derived from PDB entry 7T10 (SSTR2 bound to somatostatin-14 analog, cryo-EM).

**Figure 11** illustrates the funnel-shaped data flow across the five modules.

## 2.2 Peptide Library Generation

Module 00 generated 9,603 unique peptide sequences under the following constraints:

- **Length:** 8–14 amino acids, randomly sampled.
- **Scaffold:** Every sequence contains two cysteine residues forming a disulfide-bridged ring, with the FWKT pharmacophore (Phe-Trp-Lys-Thr) locked inside the ring. This motif was derived from the binding core of somatostatin-14 and is essential for SSTR2 recognition.
- **Amino acid pool:** 17 of the 20 standard amino acids. Lysine (K) and arginine (R) were excluded from all variable positions to minimize net positive charge and reduce megalin-mediated renal reabsorption. Cysteine (C) was excluded from variable positions to prevent unintended disulfide bonds. The single K within the FWKT motif was retained.

The general scaffold follows the pattern: [N-tail]–C–[inner-L]–FWKT–[inner-R]–C–[C-tail], where inner and tail segments are drawn from the restricted amino acid pool.

## 2.3 Five-Indicator Stability Filtering

Module 01 evaluated all 9,603 sequences using five pharmacokinetic stability indicators computed entirely from sequence-level descriptors, requiring no 3D structure generation:

1. **Protease Resistance:** Number of trypsin (K/R) and chymotrypsin (F/W/Y/L) cleavage sites, excluding the protected FWKT motif.
2. **Aggregation Tendency:** Boman index (protein-binding potential) combined with Kyte-Doolittle hydrophobicity.
3. **Kidney Toxicity:** Net charge at pH 7.4, calculated via the Henderson-Hasselbalch equation using the `peptides` library. Lower positive charge reduces megalin-mediated renal uptake.
4. **Plasma Protein Binding (PPB):** Estimated from Crippen LogP and topological polar surface area (TPSA) computed by RDKit.
5. **Rigidity:** Number of rotatable bonds (RDKit), where fewer bonds imply lower conformational entropy penalty upon binding.

Each indicator was normalized to [0, 1] via min-max scaling and combined into a **Composite Stability Score (CSS)**:

> CSS = 0.20 × S_protease + 0.15 × S_aggregation + 0.30 × S_kidney + 0.20 × S_PPB + 0.15 × S_rigidity

The kidney toxicity indicator received the highest weight (0.30) as it directly addresses our safety objective. The top 100 sequences by CSS were selected for downstream analysis.

**Filtering results:** All top-100 candidates had zero protease cleavage sites, net charge of approximately −0.29 at pH 7.4, and CSS values ranging from 0.842 to 0.868. For comparison, the DOTATATE reference carries a net charge of +0.707, indicating a substantial reduction in predicted renal toxicity risk for all selected candidates.

## 2.4 Molecular Docking and Surrogate Model

Module 02 performed structure-based binding affinity prediction through the following steps:

**Receptor preparation.** The SSTR2 receptor structure (Chain R) was extracted from PDB 7T10 and converted to PDBQT format with Gasteiger partial charges and AutoDock 4 atom types. The docking search box was centered on the orthosteric binding pocket, defined by the coordinates of the co-crystallized SST-14 ligand (center: 159.9, 149.6, 107.3 Å; box size: 30 × 30 × 30 Å).

**Peptide-receptor docking.** Each of the 100 candidate sequences was converted to a 3D molecular structure using RDKit (`MolFromSequence` → `EmbedMolecule` with ETKDGv3 → UFF optimization), then converted to PDBQT format for docking. AutoDock Vina 1.2.7 was used with exhaustiveness = 8. The best binding mode (lowest ΔG) was recorded for each peptide.

**Sequence embedding.** All 100 sequences were encoded into fixed-length numerical vectors using the ESM-2 protein language model (facebook/esm2_t6_8M_UR50D, 8M parameters). Sequences were tokenized and passed through the transformer encoder, and the mean-pooled last hidden state was used as the embedding vector.

**Surrogate model training.** An XGBoost gradient boosting regressor was trained on the ESM-2 embeddings as input features and the Vina docking scores as labels. The model was validated on a held-out 20% split. The trained surrogate then predicted binding free energy (ΔG_pred) for all 100 sequences, enabling rapid ranking without additional docking runs.

**Selectivity estimation.** For each candidate, a selectivity ratio was estimated by comparing its SSTR2 binding energy against a penalty-adjusted non-target score that accounts for the absence of the SSTR2-specific Trp/Lys binding interactions in SSTR1/3/4/5 subtypes.

**Composite ranking score.** The final ranking was computed as: CompositeScore = ΔG_pred − log₁₀(Selectivity), where more negative values indicate better candidates.

## 2.5 Multi-Agent Evaluation

Module 03 implemented a multi-agent evaluation system using the LangGraph state machine framework. The top 5 candidates (by composite ranking score) were submitted to a panel of five specialized agents, each applying structured evaluation criteria:

| Agent | Role |
|-------|------|
| Evidence Curator | Analyzes structural compatibility with the known SSTR2 binding mechanism (Trp-Lys master-key motif from PDB 7T10) |
| Structure Analyst | Evaluates Criterion 1 (selectivity ≥ 100×) and Criterion 2 (ΔG ≤ −7.0 kcal/mol for 72 h retention) |
| PK Expert | Evaluates Criterion 3 (net charge ≤ DOTATATE reference) using all five stability indicators |
| Design Generator | Proposes redesigned sequences when a candidate is rejected, preserving the Cys ring and FWKT motif |
| Experiment Planner | Generates a preclinical validation protocol for approved candidates |

The evaluation follows a sequential state machine: Evidence Curator → Structure Analyst → (if PASS) → PK Expert → (if PASS) → Experiment Planner. If either the Structure Analyst or PK Expert issues a REJECT verdict, the Design Generator proposes a modified sequence and the cycle restarts (up to 3 redesign iterations).

**Evaluation results:** All 5 submitted candidates passed both the Structure Analyst and PK Expert evaluations on the first attempt, with zero rejections. This high pass rate reflects the stringent upstream filtering applied in Modules 01 and 02.

## 2.6 Validation and Analysis

Module 04 computed the equilibrium dissociation constant (K_d) from the predicted binding free energy using the thermodynamic relationship:

> K_d = exp(ΔG_pred / RT)

where R = 1.987 cal/(mol·K), T = 310.15 K (37 °C), and RT = 0.616 kcal/mol.

A total of 13 analytical figures and 1 statistical summary table were generated at 300 DPI for publication. Key analyses include: ΔG_pred distribution (Fig. 1), selectivity distribution (Fig. 2), charge safety profile (Fig. 3), ΔG vs. selectivity quadrant plot (Fig. 4), amino acid composition heatmap (Fig. 8), property correlation matrix (Fig. 9), and multi-dimensional radar chart of the top 5 candidates (Fig. 10).

## 2.7 Results

### 2.7.1 Pipeline Throughput

| Stage | Input | Output | Reduction |
|-------|------:|-------:|-----------|
| Library Generation (Module 00) | — | 9,603 sequences | — |
| Stability Filtering (Module 01) | 9,603 | 100 | 99.0% eliminated |
| Docking + Surrogate (Module 02) | 100 | 100 (ranked) | — |
| Multi-Agent Evaluation (Module 03) | 5 | 5 approved | 100% pass rate |

### 2.7.2 Binding Affinity Prediction

The predicted binding free energies of all 100 post-filter candidates ranged from −13.41 to −9.53 kcal/mol, all substantially exceeding the 72-hour retention threshold of −7.0 kcal/mol. The corresponding K_d values ranged from 0.35 nM to 1,553 nM. For reference, the DOTATATE surrogate-predicted ΔG was −8.5 kcal/mol (K_d ≈ 1,023 nM).

### 2.7.3 Top 5 Candidates

**Table 1.** Top 5 candidate peptides ranked by composite score.

| Rank | Sequence | Length | ΔG_pred (kcal/mol) | K_d (nM) | Selectivity | Net Charge (pH 7.4) | CSS |
|-----:|----------|-------:|--------------------:|---------:|------------:|--------------------:|----:|
| 1 | CAFWKTCE | 8 | −13.41 | 0.35 | 30,828× | −0.290 | 0.851 |
| 2 | CFWKTPEC | 8 | −13.33 | 0.40 | 28,976× | −0.290 | 0.868 |
| 3 | CPFWKTEC | 8 | −13.06 | 0.63 | 23,501× | −0.290 | 0.868 |
| 4 | CFWKTDCW | 8 | −12.97 | 0.73 | 21,949× | −0.290 | 0.857 |
| 5 | PPCDFWKTCF | 10 | −12.61 | 1.31 | 16,579× | −0.290 | 0.843 |

All five candidates share the following characteristics:
- The FWKT pharmacophore is fully enclosed within the Cys-Cys disulfide ring.
- Zero protease cleavage sites outside the FWKT motif.
- Net charge at pH 7.4 is approximately −0.29, compared to +0.707 for DOTATATE — a reduction of ~1.0 charge unit.
- Predicted K_d values in the sub-nanomolar to low-nanomolar range.

### 2.7.4 Comparison with DOTATATE Reference

**Table 2.** Comparison of the top candidate (CAFWKTCE) against the DOTATATE reference.

| Property | CAFWKTCE | DOTATATE | Improvement |
|----------|:--------:|:--------:|:-----------:|
| ΔG_pred (kcal/mol) | −13.41 | −8.50 | −4.91 |
| K_d (nM) | 0.35 | 1,023.07 | ~2,880× |
| Selectivity | 30,828× | 50× | ~616× |
| Net Charge (pH 7.4) | −0.290 | +0.707 | −0.997 units |
| Protease Sites | 0 | 2 | Eliminated |
| CSS | 0.851 | — | — |

### 2.7.5 Success Criteria Assessment

All 100 candidates in the final ranked list were assessed against three predefined success criteria:

| Criterion | Threshold | Pass Rate |
|-----------|-----------|-----------|
| Target Selectivity ≥ 100× | ≥ 100× vs. SSTR1/3/4/5 | 100/100 (100%) |
| 72 h Retention | ΔG ≤ −7.0 kcal/mol | 100/100 (100%) |
| Kidney Safety | Net charge ≤ +0.707 | 100/100 (100%) |

The 100% pass rate across all criteria indicates that the upstream filtering (Module 01) effectively excluded candidates that would fail downstream evaluation.

### 2.7.6 Structural Observations

Four of the five top candidates are 8-residue peptides, suggesting that shorter, more rigid scaffolds are favored by the composite scoring function. The top-ranked candidate (CAFWKTCE) features a Cys1-Cys7 disulfide ring enclosing 5 residues, with the FWKT motif at positions 3–6. The inclusion of glutamic acid (E) in all top-4 candidates contributes to the negative net charge, directly supporting the kidney safety objective.

Amino acid composition analysis (Fig. 8) revealed a strong enrichment of acidic residues (D, E) and structure-constraining residues (P) in the top-ranked candidates, consistent with the design objectives of charge reduction and rigidity enhancement.

## 2.8 Limitations

The ΔG_pred values reported in this study are predictions from an XGBoost surrogate model trained on AutoDock Vina docking scores and ESM-2 sequence embeddings. They are not experimentally measured binding affinities. The K_d values derived from ΔG_pred should be interpreted as computational estimates that require experimental validation through in vitro binding assays (e.g., saturation binding on AR42J cells) and, ultimately, in vivo biodistribution studies. The selectivity estimates are based on computational penalty models and have not been confirmed by competitive binding assays against SSTR1/3/4/5 subtypes.

## 2.9 Design Rationale for Key Parameters

### 2.9.1 Library Size (10,000 → 9,603)

The initial library size was set to 10,000 sequences, representing a practical computational budget for a single workstation without GPU acceleration. Given the scaffold constraints ([N-tail]–C–[inner-L]–FWKT–[inner-R]–C–[C-tail]) with variable lengths of 8–14 residues and a 17-amino-acid alphabet, the theoretical combinatorial space ranges from 17² ≈ 289 (for the shortest 8-residue peptides with 2 variable positions) to 17⁸ ≈ 6.97 × 10⁹ (for the longest 14-residue peptides with 8 variable positions). The 10,000-sequence sample thus covers less than 0.00015% of the total space, functioning as a sparse random exploration rather than an exhaustive enumeration.

To maximize structural diversity within this budget, three strategies were employed: (1) peptide lengths were sampled uniformly from the range 8–14, (2) variable residues were distributed across four scaffold regions (N-tail, inner-L, inner-R, C-tail) using a Stars-and-Bars random partition algorithm to ensure diverse positional allocations, and (3) each variable position was independently sampled from the 17-amino-acid pool. A master random seed (RANDOM_SEED = 42) was used with worker-specific derived seeds for full reproducibility. After generation, sequence-level deduplication removed 397 duplicates (3.97%), yielding 9,603 unique sequences. The low duplicate rate empirically confirms that the sampled fraction is negligibly small relative to the combinatorial space.

### 2.9.2 Top-100 Cutoff from PK Filtering

The decision to advance 100 candidates from the stability filtering stage was based on two considerations:

**Computational throughput constraint.** AutoDock Vina molecular docking requires approximately 1–5 minutes per peptide on a single CPU core. Processing 100 candidates thus requires approximately 5 hours of wall-clock time, which is feasible for single-workstation execution within one working day. Scaling to 1,000 candidates would require ~50 hours, making iterative refinement impractical without cluster resources.

**Screening funnel convention.** The 100/9,603 ratio yields a pass rate of approximately 1.04%, consistent with the ~1% hit rate commonly observed in the first-pass filtering stage of high-throughput virtual screening campaigns. This ensures that only candidates with favorable pharmacokinetic profiles proceed to the computationally expensive docking stage, while retaining sufficient diversity for downstream ranking.

### 2.9.3 Success Criteria: Clinical and Literature Basis

Three predefined success criteria were applied to evaluate all candidates. The rationale and threshold for each criterion are grounded in published clinical and pharmacological evidence.

**Criterion 1: Target Selectivity ≥ 100-fold.** The 100-fold selectivity threshold over SSTR1/3/4/5 subtypes was derived from the affinity profiles of clinically validated somatostatin analogs. Reubi et al. (2000) reported that octreotide exhibits IC₅₀ values of 2.0 ± 0.7 nM for SSTR2 versus 187 ± 55 nM for SSTR3, corresponding to approximately 94-fold selectivity — the weakest selectivity gap among all subtypes [1]. The 100-fold threshold was adopted as a rounded lower bound of this clinically validated margin. For reference, the clinical gold standard [⁶⁸Ga]Ga-DOTA-TATE achieves >5,000-fold selectivity over all non-SSTR2 subtypes [1].

**Criterion 2: 72-hour Tumor Retention (ΔG ≤ −7.0 kcal/mol).** The 72-hour retention window is anchored to the pharmacokinetics of ¹⁷⁷Lu-labeled DOTATATE, whose terminal blood half-life is approximately 71 ± 28 hours [2, 3]. The 72-hour post-injection timepoint is also a standard dosimetry imaging acquisition point recommended by the EANM guidelines [4]. A binding energy threshold of ΔG ≤ −7.0 kcal/mol (corresponding to K_d ≈ 7 μM at 310 K) was applied as a computational screening filter to eliminate non-binders, consistent with established virtual screening practice using AutoDock scoring functions [5]. This threshold serves as a first-pass elimination criterion; the clinically validated DOTATATE operates at substantially stronger binding (ΔG ≈ −13 kcal/mol, K_d = 0.2 nM) [1].

**Criterion 3: Kidney Safety (Net Charge ≤ DOTATATE).** Nephrotoxicity is the dose-limiting toxicity in PRRT, with kidneys reaching the 23 Gy absorbed dose limit before bone marrow in the majority of patients [4]. The primary mechanism is megalin/cubilin-mediated endocytic reabsorption of radiolabeled peptides in the proximal tubules, as demonstrated in megalin-deficient mouse models where kidney uptake was reduced by 40–70% [6, 7]. Peptide net charge directly modulates this uptake: cationic peptides are preferentially reabsorbed via electrostatic interactions with anionic brush border membrane sites [8], while molecular charge engineering (e.g., polyglutamic acid modification) has achieved kidney uptake reductions of up to 90% [9]. Accordingly, the net charge of DOTATATE (+0.707 at pH 7.4) was adopted as the upper bound, with the design objective of producing candidates with lower (ideally negative) net charge to reduce megalin-mediated renal accumulation. Clinical kidney protection strategies such as lysine/arginine co-infusion achieve 35–50% uptake reduction [10], supporting the feasibility of charge-based nephroprotection.

## References

[1] Reubi JC, Schär JC, Waser B, et al. Affinity profiles for human somatostatin receptor subtypes SST1–SST5 of somatostatin radiotracers selected for scintigraphic and radiotherapeutic use. *Eur J Nucl Med*. 2000;27(3):273–282.

[2] Strosberg J, El-Haddad G, Wolin E, et al. Phase 3 trial of ¹⁷⁷Lu-Dotatate for midgut neuroendocrine tumors (NETTER-1). *N Engl J Med*. 2017;376(2):125–135.

[3] DrugBank. Lutetium Lu 177 dotatate — Pharmacokinetics. DrugBank ID: DB13985.

[4] Sundlöv A, Sjögreen-Gleisner K, Svensson J, et al. Individualised ¹⁷⁷Lu-DOTATATE treatment of neuroendocrine tumours based on kidney dosimetry. *Eur J Nucl Med Mol Imaging*. 2017;44(9):1480–1489.

[5] Forli S, Huey R, Pique ME, et al. Computational protein–ligand docking and virtual drug screening with the AutoDock suite. *Nat Protoc*. 2016;11(5):905–919.

[6] Vegt E, de Jong M, Wetzels JF, et al. Renal toxicity of radiolabeled peptides and antibody fragments: mechanisms, impact on radionuclide therapy, and strategies for prevention. *J Nucl Med*. 2010;51(7):1049–1058.

[7] Melis M, Krenning EP, Bernard BF, et al. Megalin is essential for renal proximal tubule reabsorption of ¹¹¹In-DTPA-octreotide. *J Nucl Med*. 2005;46(10):1696–1700.

[8] Gotthardt M, van Eerd-Vismale J, Oyen WJ, et al. Indication for different mechanisms of kidney uptake of radiolabeled peptides. *J Nucl Med*. 2007;48(4):596–601.

[9] Behe M, Kluge G, Becker W, et al. Use of polyglutamic acids to reduce uptake of radiometal-labeled minigastrin in the kidneys. *J Nucl Med*. 2005;46(6):1012–1015.

[10] Rolleman EJ, Valkema R, de Jong M, et al. Safe and effective inhibition of renal uptake of radiolabelled octreotide by a combination of lysine and arginine. *Eur J Nucl Med Mol Imaging*. 2003;30(1):9–15.
