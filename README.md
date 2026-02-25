# SSTR2 AI Scientist

**AI-Driven Discovery Pipeline for Next-Generation SSTR2-Targeted Radiopharmaceuticals (PRRT)**

> Somatostatin Receptor Type 2 (SSTR2) 표적 차세대 방사성의약품(PRRT) 후보 펩타이드를
> 대규모 생성 → 물리화학적 필터링 → AI 대리 모델 도킹 → 다중 에이전트 합의 심사까지
> 엔드투엔드로 수행하는 자동화 파이프라인

---

## Table of Contents

- [Project Overview](#project-overview)
- [3 Key Success Criteria](#3-key-success-criteria)
- [Pipeline Architecture](#pipeline-architecture)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Module Details](#module-details)
- [Results Summary](#results-summary)
- [Figures & Visualizations](#figures--visualizations)
- [Dashboard](#dashboard)
- [Validation & Limitations](#validation--limitations)
- [Configuration](#configuration)
- [License](#license)

---

## Project Overview

SSTR2(소마토스타틴 수용체 2형)는 신경내분비종양(NET) 세포 표면에 정상 세포 대비 수십 배 과발현되는 GPCR이다. 이 수용체를 정밀 타격하는 **PRRT(Peptide Receptor Radionuclide Therapy)**는 종양에 방사선을 집중 전달하여 암세포를 선택적으로 파괴하는 차세대 치료법이다.

본 프로젝트는 기존 Dotatate(Lutathera)를 능가하는 차세대 SSTR2 표적 펩타이드를 **AI 기반 파이프라인**으로 발굴한다.

### Core Innovation

```
10,000+ 펩타이드 생성 → 5대 안정성 필터 → 3D 도킹 대리 모델 → 5인 AI 드림팀 심사 → 논문용 검증
```

- **PeptiVerse/ESM-2** 단백질 언어 모델을 활용한 펩타이드 임베딩
- **XGBoost** 대리 모델로 결합 자유 에너지(ΔG) 예측
- **LangGraph** 기반 5인 AI 에이전트 합의 시스템 (Ollama 로컬 LLM)
- 학술지 규격 300 DPI 시각화 + 자동 리포트 생성

### Reference Compound: Dotatate (Lutathera)

| Property | Value |
|----------|-------|
| Sequence | D-Phe-Cys-Tyr-D-Trp-Lys-Thr-Cys-Thr (L-AA approx: `FCYWKTCT`) |
| ΔG | -8.5 kcal/mol |
| K_d | ~1,023 nM |
| Selectivity | ~50x |
| Net Charge (pH 7.4) | +0.707 |

---

## 3 Key Success Criteria

본 파이프라인의 모든 설계 결정은 아래 3대 성공 지표를 중심으로 최적화된다.

| # | Criterion | Target | Mechanism |
|---|-----------|--------|-----------|
| 1 | **100x Selectivity** | Selectivity >= 100x | FWKT pharmacophore가 SSTR2 결합 포켓에 선택적 결합. Trp이 TM3-TM7 소수성 포켓에 삽입, Lys이 Asp3.32와 염교 형성 |
| 2 | **72h Retention** | ΔG <= -7.0 kcal/mol | Cys-Cys 이황화 결합 고리 구조가 엔트로피 페널티를 줄이고, 베타-턴 구조를 안정화하여 k_off 최소화 |
| 3 | **20% Kidney Safety** | Net Charge <= +0.707 | K/R 잔기 제거로 양전하 최소화 → 신장 메갈린(Megalin) 매개 재흡수 억제 |

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                    SSTR2 AI Scientist Pipeline                      │
├─────────────┬───────────────┬───────────────┬───────────────────────┤
│  Module 00  │   Module 01   │   Module 02   │      Module 03       │
│  Generation │   PK Filter   │   Surrogate   │    AI Boardroom      │
├─────────────┼───────────────┼───────────────┼───────────────────────┤
│             │ 5-Indicator   │ PeptiVerse/   │ 5 AI Agents:         │
│ FWKT+Cys-  │ Composite     │ ESM-2 Embed   │  Evidence Curator    │
│ Cys Scaffold│ Stability     │ + XGBoost     │  Structure Analyst   │
│ 17 AA Pool │ Score         │ Surrogate     │  PK Expert           │
│ Length 8-14 │               │ ΔG Prediction │  Design Generator    │
│             │               │ + Selectivity │  Experiment Planner  │
├─────────────┼───────────────┼───────────────┼───────────────────────┤
│ 10,000 seq  │   → Top 100   │   → Top 100   │   → 4-5 Approved    │
│  (or 100M)  │  (CSS ranked) │ (Composite)   │  (Consensus Vote)   │
└─────────────┴───────────────┴───────────────┴───────────────────────┘
                                                        │
                                                        ▼
                                              ┌───────────────────┐
                                              │    Module 04      │
                                              │  Validation &     │
                                              │  Analysis         │
                                              │  (13 Figures +    │
                                              │   K_d Derivation  │
                                              │   + Report)       │
                                              └───────────────────┘
```

---

## Project Structure

```
SSTR2_AI_SCIENTIST/
│
├── config.py                  # Central configuration (all parameters)
├── 00_peptide_gen.py          # Module 00: Combinatorial peptide generation
├── 01_fast_pk_filter.py       # Module 01: 5-indicator PK stability filter
├── 02_train_surrogate.py      # Module 02: 3D docking + XGBoost surrogate
├── 03_ai_boardroom.py         # Module 03: 5-agent AI boardroom (LangGraph)
├── 04_validate_and_analyze.py # Module 04: Validation & publication figures
├── dashboard.py               # Streamlit real-time monitoring dashboard
│
├── setup_env.py               # One-click environment setup (venv + all deps)
├── setup_env.bat              # Windows launcher for setup_env.py
├── run_pipeline.py            # Sequential pipeline runner (--test-mode)
├── run_pipeline.bat           # Windows launcher for run_pipeline.py
├── launch_mission.sh          # Tmux 3-pane launcher (Linux/Mac)
├── requirements.txt           # Python dependencies (complete)
│
├── output/
│   ├── SSTR2_Raw_Library.csv              # Module 00 output (10K sequences)
│   ├── SSTR2_Filtered_Top100.csv          # Module 01 output (Top 100)
│   ├── SSTR2_Top100_Boardroom_Ready.csv   # Module 02 output (scored Top 100)
│   ├── Boardroom_Transcript.log           # Module 03 discussion transcript
│   ├── Final_Winning_Peptides_Report.md   # Module 03 final report
│   ├── Top_5_Candidates_Summary.csv       # Module 04 top-5 with K_d
│   ├── Validation_Report.md               # Module 04 validation report
│   └── figures/                           # 13 figures + tables (PNG + PDF)
│
└── Prompt_1~5.txt             # Design specification documents
```

---

## Installation

### Prerequisites

- **Python** >= 3.10
- **Ollama** (for Module 03 AI Boardroom)
- **AutoDock Vina** (optional, for production docking)

### One-Click Setup (Recommended)

```bash
# Windows
setup_env.bat

# Cross-platform
python setup_env.py
```

This automatically:
1. Creates `.venv` virtual environment
2. Detects NVIDIA GPU and installs torch (CUDA or CPU)
3. Installs all pip dependencies (requirements.txt)
4. Pre-downloads PeptiVerse & ESM-2 model weights
5. Checks Ollama installation
6. Runs full environment verification

### Setup Options

```bash
python setup_env.py                # Default (venv + GPU auto-detect + model download)
python setup_env.py --cpu-only     # Force CPU-only torch
python setup_env.py --no-venv      # Install directly (no virtual environment)
python setup_env.py --skip-models  # Skip PLM model pre-download
python setup_env.py --full         # Everything + Ollama model pull
```

### Manual Installation

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
# source .venv/bin/activate     # Linux/Mac

pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install -r requirements.txt
```

### Ollama Setup (for Module 03)

```bash
# Install: https://ollama.com/download
ollama pull llama3.1:8b
ollama serve                    # Start server (keep running)
```

---

## Quick Start

### Test Mode (No external docking engine required)

```bash
# Activate environment
.venv\Scripts\activate          # Windows
# source .venv/bin/activate     # Linux/Mac

# Start Ollama (separate terminal)
ollama serve

# Run entire pipeline
python run_pipeline.py --test-mode
```

### Pipeline Runner Options

```bash
python run_pipeline.py --test-mode                 # Test mode (dummy docking)
python run_pipeline.py --test-mode --from 02       # Restart from Module 02
python run_pipeline.py --test-mode --skip-boardroom # Skip Module 03
python run_pipeline.py --test-mode --dashboard     # Auto-launch dashboard after
```

### Individual Module Execution

```bash
python 00_peptide_gen.py                    # Generate 10K peptides
python 01_fast_pk_filter.py                 # Filter to Top 100
python 02_train_surrogate.py --test-mode    # Surrogate model + scoring
python 03_ai_boardroom.py                   # AI boardroom evaluation
python 04_validate_and_analyze.py           # Validation + figures
```

### Scale Up

```bash
NUM_SEQUENCES=100000 python 00_peptide_gen.py       # 100K sequences
NUM_SEQUENCES=100000000 python 00_peptide_gen.py     # 100M sequences
```

### Dashboard

```bash
streamlit run dashboard.py
# Open http://localhost:8501
```

---

## Module Details

### Module 00: Peptide Generation (`00_peptide_gen.py`)

대규모 SSTR2 표적 펩타이드 서열 라이브러리를 생물학적 제약 조건 하에 생성한다.

**Scaffold Architecture:**
```
[N-tail] ─ C ─ [Inner-L] ─ F─W─K─T ─ [Inner-R] ─ C ─ [C-tail]
╰─선형─╯   ╰══════ Disulfide Bond Ring (Cyclic) ══════╯   ╰─선형─╯
```

- **FWKT Pharmacophore**: SST-14 유래 결합 활성 중심부 (PDB 7T10, 7T11)
- **AA Pool**: 17 types (standard 20 minus K, R, C)
- **Length**: 8-14 residues
- **Parallelization**: `multiprocessing.Pool` for 100M-scale generation

### Module 01: PK Stability Filter (`01_fast_pk_filter.py`)

5대 물리화학적 안정성 지표를 계산하여 Composite Stability Score(CSS)로 순위를 매기고 Top-100을 추출한다.

| # | Indicator | Library | Weight | Interpretation |
|---|-----------|---------|--------|----------------|
| 1 | Protease Resistance | `re` (regex) | 0.20 | Trypsin/Chymotrypsin cleavage sites (fewer = better) |
| 2 | Aggregation Tendency | `peptides` (Boman + KD) | 0.15 | Protein interaction potential (lower = better) |
| 3 | Kidney Toxicity | `peptides` (charge pH 7.4) | 0.30 | Net charge vs Dotatate (lower = safer) |
| 4 | PPB | `RDKit` (LogP, TPSA) | 0.20 | Plasma protein binding (lower = more free drug) |
| 5 | Rigidity | `RDKit` (rotatable bonds) | 0.15 | Conformational entropy penalty (fewer = better) |

**Composite Stability Score:** CSS = w1*S_protease + w2*S_aggregation + w3*S_kidney + w4*S_ppb + w5*S_rigidity

### Module 02: Surrogate Model (`02_train_surrogate.py`)

1D 서열을 3D 구조로 변환하고, AI 대리 모델로 결합 에너지를 예측한다.

**3-Tier Embedding Cascade:**
```
1st: PeptiVerse (ChatterjeeLab/PeptiVerse)  ← Peptide-specialized PLM
2nd: ESM-2 (facebook/esm2_t6_8M_UR50D)      ← General protein PLM
3rd: Manual Feature Engineering               ← RDKit + peptides library
```

**Docking Pipeline:**
- **Production Mode**: RDKit `MolFromSequence` → `EmbedMolecule` (ETKDG) → AutoDock Vina
- **Test Mode** (`--test-mode`): SAR rule-based dummy scores + `TEST MODE` watermark

**Academic Integrity (Hard Stop Policy):**
- Production mode without Vina installed → `RuntimeError` (no fake data)
- Test mode output always carries `"TEST MODE: DUMMY DATA — NOT FOR PUBLICATION"` watermark

### Module 03: AI Boardroom (`03_ai_boardroom.py`)

LangGraph StateGraph 위에서 5인 AI 전문가가 실시간 협업하여 후보를 심사한다.

| Agent | Role | Evaluates |
|-------|------|-----------|
| Evidence Curator | Literature/PDB 기반 지식 제공 | FWKT binding mechanism validation |
| Structure Analyst | 3D 구조 심사 | 100x selectivity, 72h retention (k_off), DOTA steric hindrance |
| PK Expert | 약물동태 심사 | 5-indicator CSS, kidney filtration, PPB |
| Design Generator | 서열 재설계 | Mutations if REJECT (preserves Cys+FWKT scaffold) |
| Experiment Planner | 실험 계획서 | DOTA coupling, Lu-177 labeling, in vivo protocol |

**State Machine:**
```
Evidence Curator → Structure Analyst → [PASS] → PK Expert → [PASS] → Experiment Planner
                                        [REJECT]             [REJECT]
                                           ↓                    ↓
                                     Design Generator ←─────────┘
                                           ↓
                                     Structure Analyst (re-evaluation)
```

### Module 04: Validation & Analysis (`04_validate_and_analyze.py`)

학술지 규격의 시각화와 통계 검증을 수행한다.

- **K_d Derivation**: ΔG = RT·ln(K_d), where R = 1.987 cal/(mol·K), T = 310.15 K
- **13 Publication Figures** (300 DPI PNG + PDF vector)
- **Spearman Correlation** analysis
- **Test Mode Watermark** detection and auto-insertion

---

## Results Summary

### Pipeline Performance (Test Mode: 10K)

| Stage | Input | Output | Time |
|-------|-------|--------|------|
| Module 00: Generation | — | 9,602 unique sequences | 1.40s |
| Module 01: PK Filter | 9,602 | 100 (CSS Top-100) | 3.38s |
| Module 02: Surrogate | 100 | 100 (scored + ranked) | ~10s |
| Module 03: Boardroom | 5 evaluated | 4 approved | ~2min |
| Module 04: Analysis | 100 | 13 figures + report | ~5s |

### Top 5 Candidates

| Rank | Sequence | ΔG_pred (kcal/mol) | K_d (nM) | Selectivity | Net Charge |
|------|----------|--------------------|----------|-------------|------------|
| 1 | `PCDPFWKTCW` | -8.244 | 1,549.9 | 1,166.3x | -0.290 |
| 2 | `ECFWKTICP` | -8.169 | 1,750.5 | 1,310.4x | -0.290 |
| 3 | `ACPFWKTCD` | -8.062 | 2,082.4 | 633.4x | -0.290 |
| 4 | `CDFWKTPC` | -8.305 | 1,403.9 | 360.8x | -0.290 |
| 5 | `CFWKTPEC` | -8.056 | 2,102.8 | 499.6x | -0.290 |

### Success Criteria Achievement (Top 100)

| Criterion | Threshold | Pass Rate |
|-----------|-----------|-----------|
| 100x Selectivity | >= 100x | **97/100 (97.0%)** |
| 72h Retention | <= -7.0 kcal/mol | **83/100 (83.0%)** |
| Kidney Safety | <= +0.707 charge | **100/100 (100.0%)** |

### vs Dotatate Improvement

- **Net Charge**: average -0.290 vs Dotatate +0.707 → **140% improvement**
- **Selectivity**: average ~800x vs Dotatate ~50x → **16x improvement**
- **All candidates**: 0 protease cleavage sites (Dotatate has multiple)

---

## Figures & Visualizations

All figures are saved in `output/figures/` at 300 DPI (PNG) with PDF vector copies.

| Figure | Description | Success Criteria |
|--------|-------------|-----------------|
| Fig 01 | ΔG_pred Distribution (Histogram + KDE) | Criterion 2 |
| Fig 02 | Selectivity Ratio Analysis (Box + Strip) | Criterion 1 |
| Fig 03 | Net Charge Safety Profile vs Dotatate | Criterion 3 |
| Fig 04 | ΔG vs Selectivity Quadrant (Winning Zone) | Criterion 1 + 2 |
| Fig 05 | Physicochemical Property Space (Hydro vs pI) | — |
| Fig 06 | Peptide Length vs ΔG (Box Plot) | — |
| Fig 07 | Composite Stability Score Ranking (Top-20) | — |
| Fig 08 | Amino Acid Positional Frequency Heatmap | — |
| Fig 09 | Property Correlation Matrix | — |
| Fig 10 | Radar Chart: Top-5 vs Dotatate | All criteria |
| Fig 11 | Pipeline Funnel (Sequence Compression) | — |
| Fig 12 | Spearman Correlation (Surrogate Validation) | — |
| Fig 13 | 5-Indicator Stability Profile (Top-10) | — |
| Table 01 | Descriptive Statistics Summary | — |

---

## Dashboard

Streamlit 기반 실시간 Mission Control 대시보드.

```bash
streamlit run dashboard.py
# http://localhost:8501
```

**Features:**
- **Top**: Pipeline 4-stage progress bar (file existence detection)
- **Sidebar**: Agent Monitor (real-time Transcript.log parsing with PASS/REJECT badges)
- **Center**: Candidate Ranking table with 3-criteria PASS/FAIL badges
- **Bottom**: Interactive Plotly charts (Binding Energy / Selectivity / Net Charge)
- **Auto-refresh**: Optional 5-second polling

---

## Validation & Limitations

### Terminology

All binding energies in this pipeline are reported as **AI-predicted Binding Free Energy (ΔG_pred)**, computed by an XGBoost surrogate model trained on PLM embeddings and RDKit molecular descriptors.

### K_d Derivation

```
ΔG = RT · ln(K_d)
K_d = exp(ΔG / RT)

R = 1.987 cal/(mol·K)
T = 310.15 K (37°C)
RT = 0.616 kcal/mol
```

### Disclaimer

> This Surrogate Model combines PeptiVerse/ESM-2 embeddings, RDKit molecular descriptors, and XGBoost regression. The reported ΔG_pred values are AI-predicted approximations, **NOT experimentally measured binding affinities**. All final candidates MUST be validated through physical molecular docking (e.g., AutoDock Vina, Glide) and molecular dynamics (MD) simulations before any experimental synthesis or in vivo testing.

### Future Validation Steps

1. **Physical Docking**: AutoDock Vina / Glide with SSTR2 crystal structure (PDB 7T10)
2. **MD Simulation**: OpenMM / GROMACS for binding stability (100ns+ trajectories)
3. **In Vitro**: AR42J cell binding assay, CHO-SSTR1/3/4/5 selectivity panel
4. **Radiolabeling**: DOTA-NHS coupling + Lu-177 labeling (target: >95% RCP)
5. **In Vivo**: Nude mouse AR42J xenograft biodistribution (1h-72h timepoints)

---

## Configuration

All parameters are centralized in `config.py` and overridable via environment variables.

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NUM_SEQUENCES` | 10,000 | Number of peptides to generate |
| `FILTER_TOP_N` | 100 | Top-N from PK filter |
| `PEPTIVERSE_MODEL_NAME` | ChatterjeeLab/PeptiVerse | 1st-priority PLM |
| `ESM2_MODEL_NAME` | facebook/esm2_t6_8M_UR50D | 2nd-priority PLM |
| `OLLAMA_MODEL` | llama3.1:8b | LLM for AI boardroom |
| `BOARDROOM_TOP_N` | 5 | Candidates evaluated in boardroom |
| `MAX_REDESIGN_ITERATIONS` | 3 | Max mutation cycles per candidate |
| `TRAINING_SUBSET_SIZE` | 500 | Docking training subset |
| `TOP_K_CANDIDATES` | 100 | Final output size |

### Stability Score Weights

| Weight | Value | Indicator |
|--------|-------|-----------|
| `W_PROTEASE` | 0.20 | Protease Resistance |
| `W_AGGREGATION` | 0.15 | Aggregation Tendency |
| `W_KIDNEY` | 0.30 | Kidney Toxicity (highest weight) |
| `W_PPB` | 0.20 | Plasma Protein Binding |
| `W_RIGIDITY` | 0.15 | Structural Rigidity |

### Environment Variable Override

```bash
NUM_SEQUENCES=100000000 OLLAMA_MODEL=llama3:70b python run_pipeline.py
```

---

## Tech Stack

| Layer | Technology |
|-------|------------|
| Peptide Generation | NumPy, multiprocessing |
| PK Analysis | peptides, RDKit |
| PLM Embedding | PeptiVerse, ESM-2 (HuggingFace Transformers, PyTorch) |
| Surrogate Model | XGBoost, scikit-learn |
| AI Agents | LangGraph, LangChain, Ollama |
| Visualization | Matplotlib, Seaborn, Plotly |
| Dashboard | Streamlit |
| Statistics | SciPy (Spearman correlation) |

---

## License

This project is for research and educational purposes. The AI-predicted results require experimental validation before any clinical application.
