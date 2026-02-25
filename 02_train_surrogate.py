#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  Module 02: 1D→3D 변환 + 분자 도킹 + AI 대리 모델
  (02_train_surrogate.py)
==========================================================================

  [Stage]  Phase 2 - Step 1/2
  [Input]  output/SSTR2_Filtered_Top100.csv (Module 01 출력)
  [Output] output/SSTR2_Top100_Boardroom_Ready.csv (AI 드림팀 입력)

  [학술적 무결성 원칙 — No Hallucination]
  본 모듈은 절대로 가짜 도킹 점수를 생성하지 않습니다.
  - Production Mode: AutoDock Vina 미설치 시 RuntimeError → Hard Stop
  - Test Mode (--test-mode): 명시적 더미 데이터 + 워터마크 경고

  [3D 도킹 파이프라인]
  1. 수용체(Receptor): 로컬 SSTR2_receptor.pdbqt 파일 고정
  2. 리간드(Peptide): RDKit MolFromSequence → EmbedMolecule (ETKDG) → 3D .pdb
  3. 도킹: subprocess.run으로 AutoDock Vina 호출 → ΔG 파싱
  4. 대리 모델: ESM-2/수동 임베딩 + XGBoost → 나머지 서열 ΔG 예측
  5. 선택성 + 복합 점수 → Top-100 추출

  [Usage]
  Production:  python 02_train_surrogate.py
  Test Mode:   python 02_train_surrogate.py --test-mode
==========================================================================
"""

import argparse
import numpy as np
import pandas as pd
import subprocess
import sys
import tempfile
import time
import warnings
from pathlib import Path

# -- XGBoost & scikit-learn --
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_absolute_error

# -- Peptides --
from peptides import Peptide

# -- RDKit (1D→3D 변환 필수) --
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# -- Meeko (Vina 1.2.x 공식 PDBQT 도구) --
try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    MEEKO_AVAILABLE = True
except ImportError:
    MEEKO_AVAILABLE = False

# -- ESM-2 조건부 임포트 --
try:
    import torch
    from transformers import AutoTokenizer, AutoModel
    ESM2_AVAILABLE = True
except ImportError:
    ESM2_AVAILABLE = False

# -- 프로젝트 중앙 설정 --
from config import (
    OUTPUT_DIR, DOTATATE_SEQ, NUM_WORKERS, RANDOM_SEED,
    PEPTIVERSE_MODEL_NAME, ESM2_MODEL_NAME,
    TRAINING_SUBSET_SIZE, TOP_K_CANDIDATES,
    DOCKING_RANDOM_SEED,
    RECEPTOR_PATH, VINA_EXECUTABLE,
    DOCK_CENTER_X, DOCK_CENTER_Y, DOCK_CENTER_Z,
    DOCK_SIZE_X, DOCK_SIZE_Y, DOCK_SIZE_Z,
    DOCK_EXHAUSTIVENESS,
    RT_BODY_KCAL,
)

warnings.filterwarnings("ignore", category=FutureWarning)

# Windows UTF-8
if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# Global test mode flag (set via argparse)
TEST_MODE = False

# ─── 테스트 모드 워터마크 ─────────────────────────────────────────────────
WATERMARK = "⚠ TEST MODE: DUMMY DATA — NOT FOR PUBLICATION ⚠"

# ─── 물리 상수 ────────────────────────────────────────────────────────────
DOTATATE_DOCKING_SCORE = -8.5  # kcal/mol (literature reference)
MIN_RETENTION_ENERGY = -7.0    # kcal/mol (72h retention threshold)


# ==========================================================================
#  Section 1: 1D→3D 변환 + PDBQT 준비
# ==========================================================================

def _sequence_to_mol3d(sequence):
    """
    아미노산 서열 → RDKit Mol (수소 포함, 3D 좌표 내장).

    Returns
    -------
    Chem.Mol : 3D conformer + hydrogens 포함 Mol 객체
    """
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit required. Install: pip install rdkit-pypi")

    mol = Chem.MolFromSequence(sequence)
    if mol is None:
        mol = Chem.MolFromFASTA(f">peptide\n{sequence}")
    if mol is None:
        raise ValueError(f"Cannot convert sequence to Mol: {sequence}")

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = DOCKING_RANDOM_SEED
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            raise RuntimeError(f"3D embedding failed: {sequence}")

    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass

    return mol


def _mol_to_pdbqt(mol, pdbqt_path):
    """
    RDKit Mol → PDBQT (Vina 1.2.x 호환).

    1순위: meeko (Vina 공식 PDBQT 도구, 정확한 torsion tree 생성)
    2순위: 간이 변환 (ROOT/ENDROOT/TORSDOF 0, 강체 도킹)
    """
    # ── 1순위: meeko ──
    if MEEKO_AVAILABLE:
        try:
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(
                mol_setups[0])
            if is_ok:
                with open(pdbqt_path, "w") as f:
                    f.write(pdbqt_string)
                return
        except Exception:
            pass  # Fallback to simple method

    # ── 2순위: 간이 변환 (강체 도킹) ──
    mol_no_h = Chem.RemoveHs(mol)
    pdb_block = Chem.MolToPDBBlock(mol_no_h)

    with open(pdbqt_path, "w") as f:
        f.write("ROOT\n")
        for line in pdb_block.split("\n"):
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) >= 78 else ""
                if not element:
                    for ch in atom_name:
                        if ch.isalpha():
                            element = ch.upper()
                            break
                if not element:
                    element = "C"
                # AD4 atom type
                ad4 = {"C": "C", "N": "N", "O": "OA", "S": "SA",
                       "H": "HD", "P": "P"}.get(element.upper(), "C")
                pdb_part = line[:66].ljust(66)
                f.write(f"{pdb_part}    +0.0000 {ad4:<2s}\n")
        f.write("ENDROOT\n")
        f.write("TORSDOF 0\n")


def sequence_to_3d_pdb(sequence, output_path):
    """하위 호환용 래퍼: 서열 → PDB 파일."""
    mol = _sequence_to_mol3d(sequence)
    mol_no_h = Chem.RemoveHs(mol)
    Chem.MolToPDBFile(mol_no_h, str(output_path))
    return True


# ==========================================================================
#  Section 2: AutoDock Vina 도킹 래퍼
# ==========================================================================

def _check_vina_installed():
    """AutoDock Vina 실행파일이 시스템에 존재하는지 확인한다."""
    try:
        result = subprocess.run(
            [VINA_EXECUTABLE, "--version"],
            capture_output=True, text=True, timeout=10
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def dock_peptide_vina(ligand_pdb_or_seq, receptor_pdbqt, sequence=None):
    """
    AutoDock Vina로 펩타이드-수용체 도킹을 수행하고 ΔG를 파싱한다.

    Vina 1.2.x 호환: meeko PDBQT 생성 + stdout 파싱.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        ligand_pdbqt = tmpdir / "ligand.pdbqt"
        output_pdbqt = tmpdir / "output.pdbqt"

        # 리간드 PDBQT 준비
        if sequence is not None:
            # 서열에서 직접 3D Mol → PDBQT (meeko 우선)
            mol = _sequence_to_mol3d(sequence)
            _mol_to_pdbqt(mol, ligand_pdbqt)
        else:
            # 기존 PDB 파일에서 변환 (하위 호환)
            _pdb_to_pdbqt_simple(str(ligand_pdb_or_seq), ligand_pdbqt)

        # Vina 실행
        cmd = [
            VINA_EXECUTABLE,
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", str(DOCK_CENTER_X),
            "--center_y", str(DOCK_CENTER_Y),
            "--center_z", str(DOCK_CENTER_Z),
            "--size_x", str(DOCK_SIZE_X),
            "--size_y", str(DOCK_SIZE_Y),
            "--size_z", str(DOCK_SIZE_Z),
            "--exhaustiveness", str(DOCK_EXHAUSTIVENESS),
            "--out", str(output_pdbqt),
        ]

        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300
        )

        if result.returncode != 0:
            err_detail = (result.stderr.strip() + " " +
                          result.stdout.strip()).strip()
            raise RuntimeError(
                f"Vina failed (exit code {result.returncode}): "
                f"{err_detail[:500]}"
            )

        return _parse_vina_stdout(result.stdout)


def _pdb_to_pdbqt_simple(pdb_path, pdbqt_path):
    """간이 PDB → PDBQT (ROOT/ENDROOT/TORSDOF 포함, Vina 1.2.x 호환)."""
    with open(pdb_path, "r") as f_in:
        lines = f_in.readlines()

    with open(pdbqt_path, "w") as f_out:
        f_out.write("ROOT\n")
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) >= 78 else ""
                if not element:
                    for ch in atom_name:
                        if ch.isalpha():
                            element = ch.upper()
                            break
                if not element:
                    element = "C"
                ad4 = {"C": "C", "N": "N", "O": "OA", "S": "SA",
                       "H": "HD", "P": "P"}.get(element.upper(), "C")
                pdb_part = line[:66].ljust(66)
                f_out.write(f"{pdb_part}    +0.0000 {ad4:<2s}\n")
        f_out.write("ENDROOT\n")
        f_out.write("TORSDOF 0\n")


def _parse_vina_stdout(stdout_text):
    """Vina 1.2.x stdout에서 best binding affinity (mode 1)를 파싱한다."""
    for line in stdout_text.split("\n"):
        line = line.strip()
        parts = line.split()
        if len(parts) >= 4 and parts[0] == "1":
            try:
                return float(parts[1])
            except ValueError:
                continue
    raise RuntimeError(
        f"Cannot parse Vina output. stdout={stdout_text[:300]}"
    )


# ==========================================================================
#  Section 3: 도킹 오케스트레이터 (Production vs Test Mode)
# ==========================================================================

def dock_sequences(sequences, indices):
    """
    선택된 서열들에 대해 도킹을 수행하고 ΔG 값을 반환한다.

    Production Mode: Vina 실제 도킹 수행
    Test Mode: 명시적 더미 데이터 반환 (워터마크 포함)
    """
    global TEST_MODE

    if TEST_MODE:
        print(f"\n  {WATERMARK}")
        print(f"  Test Mode: 더미 도킹 점수를 SAR 규칙 기반으로 생성합니다.")
        print(f"  이 데이터는 학술 출판에 사용할 수 없습니다.")
        rng = np.random.RandomState(DOCKING_RANDOM_SEED)
        scores = []
        for idx in indices:
            score = _generate_dummy_docking_score(sequences[idx], rng)
            scores.append(score)
        return np.array(scores)

    # ── Production Mode ──
    # 1) 수용체 파일 확인
    receptor = Path(RECEPTOR_PATH)
    if not receptor.exists():
        raise RuntimeError(
            f"학술적 무결성을 위해 수용체 파일이 필요합니다.\n"
            f"  경로: {receptor.resolve()}\n"
            f"  SSTR2_receptor.pdbqt 파일을 프로젝트 루트에 배치하세요.\n"
            f"  (PDB 7T10에서 prepare_receptor4.py로 생성 가능)\n"
            f"  테스트 실행: python 02_train_surrogate.py --test-mode"
        )

    # 2) Vina 설치 확인
    if not _check_vina_installed():
        raise RuntimeError(
            f"학술적 무결성을 위해 도킹 엔진(AutoDock Vina 등) 설치가 필수입니다.\n"
            f"가짜 데이터를 생성할 수 없습니다.\n\n"
            f"  설치 방법: https://vina.scripps.edu/\n"
            f"  또는: conda install -c conda-forge autodock-vina\n"
            f"  테스트 실행: python 02_train_surrogate.py --test-mode"
        )

    # 3) 실제 도킹 수행 (meeko + Vina 1.2.x)
    print(f"  PDBQT 준비: {'meeko (torsion tree)' if MEEKO_AVAILABLE else '간이 변환 (rigid)'}")
    scores = []
    total = len(indices)
    dock_start = time.time()
    for i, idx in enumerate(indices):
        seq = sequences[idx]
        t0 = time.time()

        try:
            dg = dock_peptide_vina(None, receptor, sequence=seq)
            scores.append(dg)
            elapsed_one = time.time() - t0
            elapsed_total = time.time() - dock_start
            avg_per_seq = elapsed_total / (i + 1)
            remaining = avg_per_seq * (total - i - 1)
            print(f"  [{i+1}/{total}] {seq:<16s} ΔG={dg:>7.3f} kcal/mol "
                  f"({elapsed_one:.1f}s, 남은 예상: {remaining/60:.1f}min)",
                  flush=True)
        except Exception as e:
            print(f"  [{i+1}/{total}] {seq:<16s} [FAIL] {e}", flush=True)
            scores.append(np.nan)

    print(f"  도킹 완료: {total}개, 총 {(time.time()-dock_start)/60:.1f}분")
    return np.array(scores)


def _generate_dummy_docking_score(sequence, rng):
    """
    TEST MODE 전용: SAR 규칙 기반 더미 도킹 점수 생성.
    이 데이터는 학술 출판에 사용할 수 없습니다.
    """
    score = -5.0
    c_pos = [i for i, aa in enumerate(sequence) if aa == "C"]
    if len(c_pos) >= 2:
        ring_size = c_pos[1] - c_pos[0] - 1
        if 5 <= ring_size <= 7:
            score -= 2.0
        elif 4 <= ring_size <= 8:
            score -= 1.0
        else:
            score += 0.5

    w_idx = sequence.find("W")
    if w_idx >= 0:
        neighbors = ""
        if w_idx > 0:
            neighbors += sequence[w_idx - 1]
        if w_idx < len(sequence) - 1:
            neighbors += sequence[w_idx + 1]
        aromatic_neighbors = sum(1 for aa in neighbors if aa in "FYW")
        score -= aromatic_neighbors * 0.4

    if len(c_pos) >= 2:
        ring_seq = sequence[c_pos[0]+1:c_pos[1]]
        acidic_in_ring = sum(1 for aa in ring_seq if aa in "DE")
        score += acidic_in_ring * 0.5
        hydro_in_ring = sum(1 for aa in ring_seq if aa in "FWYILV")
        score -= hydro_in_ring * 0.2
        if "P" in ring_seq:
            score -= 0.4

    if len(sequence) > 12:
        score += 0.25 * (len(sequence) - 12)

    score += rng.normal(0, 0.4)
    return round(score, 3)


# ==========================================================================
#  Section 4: 선택성 계산
# ==========================================================================

def calculate_selectivity_batch(sstr2_scores, sequences, rng):
    """
    비표적 수용체(SSTR1/3/4/5) 대비 선택성 비율을 계산한다.

    selectivity = Kd_nontarget / Kd_SSTR2 = exp(DDG / RT)

    [Test Mode] 비표적 점수는 SSTR2 점수의 구조적 약화 비율로 추정.
    [Production] 실제로는 비표적 수용체에도 개별 도킹을 수행해야 하나,
    현 파이프라인에서는 SAR 기반 약화 모델을 사용한다.
    이 한계점은 04_validate_and_analyze.py의 Disclaimer에서 명시.
    """
    receptors = [
        ("SSTR1", 0.55),
        ("SSTR3", 0.50),
        ("SSTR4", 0.35),
        ("SSTR5", 0.70),
    ]

    selectivities = []
    for i, seq in enumerate(sequences):
        sstr2_sc = sstr2_scores[i]
        nt_scores = []
        for name, factor in receptors:
            nt = sstr2_sc * factor
            nt = min(nt, -0.5)
            nt_scores.append(nt)

        avg_nt = np.mean(nt_scores)
        ddg = avg_nt - sstr2_sc
        sel = np.exp(ddg / RT_BODY_KCAL)
        selectivities.append(round(sel, 1))

    return np.array(selectivities)


# ==========================================================================
#  Section 5: 서열 임베딩
# ==========================================================================

def _compute_peptiverse_embeddings(sequences, model_name, batch_size=16):
    """PeptiVerse 펩타이드 특화 언어 모델 임베딩 (1순위).

    ChatterjeeLab/PeptiVerse: 펩타이드 기능 예측에 최적화된 PLM.
    Mean Pooling (attention mask 가중) 방식으로 서열 수준 임베딩을 추출한다.
    """
    print(f"  PeptiVerse 모델 로드 중: {model_name}")
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True)
    model.eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    print(f"  디바이스: {device}")

    all_embeddings = []
    total = len(sequences)
    for i in range(0, total, batch_size):
        batch = sequences[i:i+batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True,
                           truncation=True, max_length=64)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        with torch.no_grad():
            outputs = model(**inputs)
        # Mean Pooling: attention_mask로 패딩 토큰을 제외한 평균
        attention_mask = inputs["attention_mask"].unsqueeze(-1).float()
        masked = outputs.last_hidden_state * attention_mask
        mean_emb = masked.sum(dim=1) / attention_mask.sum(dim=1).clamp(min=1e-9)
        all_embeddings.append(mean_emb.cpu().numpy())
        if (i // batch_size) % 20 == 0:
            pct = min(i + batch_size, total) / total * 100
            print(f"\r  PeptiVerse 임베딩: {pct:.0f}%", end="", flush=True)
    print()
    return np.vstack(all_embeddings)


def _compute_esm2_embeddings(sequences, model_name, batch_size=16):
    """ESM-2 단백질 언어 모델 임베딩 (2순위 Fallback)."""
    print(f"  ESM-2 모델 로드 중: {model_name}")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    print(f"  디바이스: {device}")

    all_embeddings = []
    total = len(sequences)
    for i in range(0, total, batch_size):
        batch = sequences[i:i+batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True,
                           truncation=True, max_length=20)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        with torch.no_grad():
            outputs = model(**inputs)
        attention_mask = inputs["attention_mask"].unsqueeze(-1).float()
        masked = outputs.last_hidden_state * attention_mask
        mean_emb = masked.sum(dim=1) / attention_mask.sum(dim=1).clamp(min=1e-9)
        all_embeddings.append(mean_emb.cpu().numpy())
        if (i // batch_size) % 20 == 0:
            pct = min(i + batch_size, total) / total * 100
            print(f"\r  ESM-2 임베딩: {pct:.0f}%", end="", flush=True)
    print()
    return np.vstack(all_embeddings)


def _compute_manual_embeddings(sequences):
    """수동 특성 공학 기반 임베딩 (3순위 최종 Fallback)."""
    AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
    features = []
    for seq in sequences:
        feat = []
        seq_len = max(len(seq), 1)
        for aa in AA_LIST:
            feat.append(seq.count(aa) / seq_len)
        pep = Peptide(seq)
        feat.append(pep.charge(pH=7.4))
        feat.append(pep.isoelectric_point())
        feat.append(pep.hydrophobicity("Eisenberg"))
        try:
            feat.append(pep.molecular_weight())
        except Exception:
            feat.append(seq_len * 110.0)
        try:
            feat.append(pep.boman())
        except Exception:
            feat.append(0.0)
        feat.append(float(seq_len))
        c_pos = [i for i, a in enumerate(seq) if a == "C"]
        if len(c_pos) >= 2:
            ring_size = c_pos[1] - c_pos[0] - 1
            feat.append(ring_size / seq_len)
            ring_seq = seq[c_pos[0]+1:c_pos[1]]
            hydro_r = sum(1 for a in ring_seq if a in "FWYILV") / max(len(ring_seq), 1)
            feat.append(hydro_r)
            tail_r = (c_pos[0] + (seq_len - c_pos[1] - 1)) / seq_len
            feat.append(tail_r)
        else:
            feat.extend([0.0, 0.0, 0.0])
        fwkt_idx = seq.find("FWKT")
        feat.append(fwkt_idx / seq_len if fwkt_idx >= 0 else 0.5)
        features.append(feat)
    return np.array(features, dtype=np.float32)


def compute_embeddings(sequences):
    """서열 임베딩 — 3단 Cascade Fallback.

    1순위: PeptiVerse (펩타이드 특화 PLM)
    2순위: ESM-2 (범용 단백질 PLM)
    3순위: 수동 특성 공학 (torch/transformers 불필요)
    """
    if ESM2_AVAILABLE:
        # ── 1순위: PeptiVerse ──
        try:
            print("  [모드] PeptiVerse 펩타이드 특화 언어 모델 임베딩")
            return _compute_peptiverse_embeddings(
                sequences, PEPTIVERSE_MODEL_NAME)
        except Exception as e:
            print(f"  [경고] PeptiVerse 실패 ({e})")

        # ── 2순위: ESM-2 ──
        try:
            print("  [Fallback] ESM-2 단백질 언어 모델 임베딩으로 전환")
            return _compute_esm2_embeddings(sequences, ESM2_MODEL_NAME)
        except Exception as e:
            print(f"  [경고] ESM-2 실패 ({e})")

        # ── 3순위: 수동 특성 공학 ──
        print("  [Fallback] 수동 특성 공학 임베딩으로 전환")
        return _compute_manual_embeddings(sequences)
    else:
        print("  [모드] 수동 특성 공학 임베딩 (torch/transformers 미설치)")
        return _compute_manual_embeddings(sequences)


# ==========================================================================
#  Section 6: XGBoost 대리 모델
# ==========================================================================

def train_surrogate(X_train, y_train, X_val, y_val):
    """XGBoost 대리 모델 학습 및 검증."""
    model = xgb.XGBRegressor(
        n_estimators=300, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        reg_alpha=0.1, reg_lambda=1.0,
        random_state=RANDOM_SEED, n_jobs=NUM_WORKERS, verbosity=0,
    )
    model.fit(X_train, y_train, eval_set=[(X_val, y_val)], verbose=False)

    y_pred = model.predict(X_val)
    r2 = r2_score(y_val, y_pred)
    mae = mean_absolute_error(y_val, y_pred)
    print(f"  XGBoost 검증: R²={r2:.4f}, MAE={mae:.4f} kcal/mol")
    return model


# ==========================================================================
#  Main Pipeline
# ==========================================================================

def main():
    global TEST_MODE

    parser = argparse.ArgumentParser(description="SSTR2 AI Surrogate Model")
    parser.add_argument("--test-mode", action="store_true",
                        help="Run with dummy docking data (NOT for publication)")
    args = parser.parse_args()
    TEST_MODE = args.test_mode

    print("=" * 70)
    print("  SSTR2 AI Scientist — Module 02: 3D Docking + Surrogate Model")
    if TEST_MODE:
        print(f"  {WATERMARK}")
    print("=" * 70)

    start_time = time.time()

    # ── Step 1: 데이터 로드 ──
    input_path = OUTPUT_DIR / "SSTR2_Filtered_Top100.csv"

    # 하위 호환: 구 파일명도 지원
    if not input_path.exists():
        legacy_path = OUTPUT_DIR / "SSTR2_Filtered_10K.csv"
        if legacy_path.exists():
            input_path = legacy_path
        else:
            print(f"  [오류] 입력 파일 없음: {input_path}")
            sys.exit(1)

    output_path = OUTPUT_DIR / "SSTR2_Top100_Boardroom_Ready.csv"

    print(f"\n>>> Step 1: 필터 통과 서열 로드")
    df = pd.read_csv(input_path, encoding="utf-8")
    sequences = df["Sequence"].tolist()
    n_total = len(sequences)
    print(f"  로드: {n_total:,}개 서열 ({input_path.name})")

    # ── Step 2: 임베딩 ──
    print(f"\n>>> Step 2: 서열 임베딩 계산")
    embeddings = compute_embeddings(sequences)
    print(f"  임베딩: {embeddings.shape[1]}D × {n_total}개")

    # ── Step 3: 도킹 (훈련 서브셋) ──
    subset_size = min(TRAINING_SUBSET_SIZE, n_total)
    print(f"\n>>> Step 3: 분자 도킹 ({subset_size}개 서브셋)")
    if TEST_MODE:
        print(f"  [TEST MODE] SAR 규칙 기반 더미 점수 생성")
    else:
        print(f"  [PRODUCTION] AutoDock Vina 실제 도킹")
        print(f"  수용체: {RECEPTOR_PATH}")

    rng = np.random.RandomState(DOCKING_RANDOM_SEED)
    train_indices = rng.choice(n_total, size=subset_size, replace=False)

    train_scores = dock_sequences(sequences, train_indices)

    # NaN 처리 (도킹 실패 서열 제거)
    valid_mask = ~np.isnan(train_scores)
    train_indices = train_indices[valid_mask]
    train_scores = train_scores[valid_mask]

    print(f"  유효 도킹 결과: {len(train_scores)}/{subset_size}")
    print(f"  점수 범위: [{train_scores.min():.3f}, {train_scores.max():.3f}] kcal/mol")

    # ── Step 4: XGBoost 학습 ──
    print(f"\n>>> Step 4: XGBoost 대리 모델 학습")
    X_train_full = embeddings[train_indices]
    X_tr, X_va, y_tr, y_va = train_test_split(
        X_train_full, train_scores, test_size=0.2, random_state=RANDOM_SEED
    )
    print(f"  훈련: {len(X_tr)} / 검증: {len(X_va)}")
    model = train_surrogate(X_tr, y_tr, X_va, y_va)

    # ── Step 5: 전체 예측 ──
    print(f"\n>>> Step 5: 전체 서열 ΔG_pred 예측")
    predicted_scores = model.predict(embeddings)
    print(f"  예측 범위: [{predicted_scores.min():.3f}, {predicted_scores.max():.3f}]")

    # ── Step 6: 선택성 ──
    print(f"\n>>> Step 6: 선택성 계산")
    selectivities = calculate_selectivity_batch(predicted_scores, sequences, rng)
    print(f"  선택성 범위: [{selectivities.min():.1f}x, {selectivities.max():.1f}x]")

    # ── Step 7: 복합 점수 + Top-K ──
    print(f"\n>>> Step 7: 복합 점수 → Top-{TOP_K_CANDIDATES}")
    log_sel = np.log10(np.clip(selectivities, 1.0, None))
    composite = predicted_scores - log_sel  # 더 음수 = 더 좋음

    top_k = min(TOP_K_CANDIDATES, n_total)
    top_indices = np.argsort(composite)[:top_k]

    results = []
    for rank, idx in enumerate(top_indices, 1):
        row = df.iloc[idx]
        results.append({
            "Rank": rank,
            "Sequence": row["Sequence"],
            "Length": row.get("Length", len(row["Sequence"])),
            "SSTR2_BindingEnergy": round(predicted_scores[idx], 3),
            "Selectivity": round(selectivities[idx], 1),
            "CompositeScore": round(composite[idx], 3),
            "NetCharge": round(row.get("NetCharge", 0), 4),
            "pI": round(row.get("pI", 0), 4),
            "Hydrophobicity": round(row.get("Hydrophobicity", 0), 4),
            # Module 01의 5대 지표 전달 (Boardroom용)
            "ProteaseSites": int(row.get("ProteaseSites", 0)),
            "AggregationScore": round(row.get("AggregationScore", 0), 4),
            "PPBScore": round(row.get("PPBScore", 0), 4),
            "Rigidity": int(row.get("Rigidity", 0)),
            "CompositeStabilityScore": round(
                row.get("CompositeStabilityScore", 0), 4),
        })

    df_top = pd.DataFrame(results)

    # TEST MODE일 경우 WATERMARK 컬럼 추가 (04 모듈에서 감지용)
    if TEST_MODE:
        df_top["WATERMARK"] = WATERMARK

    # ── Step 8: 저장 ──
    print(f"\n>>> Step 8: 결과 저장")
    df_top.to_csv(output_path, index=False, encoding="utf-8")
    if TEST_MODE:
        print(f"  {WATERMARK}")
    print(f"  출력: {output_path}")

    elapsed = time.time() - start_time

    print(f"\n{'=' * 70}")
    print(f"  >>> 대리 모델 파이프라인 완료!")
    if TEST_MODE:
        print(f"  {WATERMARK}")
    print(f"  입력: {n_total:,} → 출력: {len(df_top):,} (Top-{top_k})")
    print(f"  소요 시간: {elapsed:.2f}초")
    print(f"{'=' * 70}")

    # Top 5
    print(f"\n  [Top 5 후보]")
    print(f"  {'Rank':>4} | {'Sequence':<16} | {'ΔG_pred':>8} | "
          f"{'Sel':>7} | {'Charge':>7}")
    print(f"  {'-'*4}-+-{'-'*16}-+-{'-'*8}-+-{'-'*7}-+-{'-'*7}")
    for _, row in df_top.head(5).iterrows():
        print(f"  {int(row['Rank']):>4} | {row['Sequence']:<16} | "
              f"{row['SSTR2_BindingEnergy']:>8.3f} | "
              f"{row['Selectivity']:>6.1f}x | "
              f"{row['NetCharge']:>+6.3f}")

    return output_path


if __name__ == "__main__":
    main()
