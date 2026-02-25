#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
  SSTR2 AI Scientist — 파이프라인 순차 실행기  (run_pipeline.py)
==========================================================================

  5개 모듈을 순차적으로 실행하고, 각 단계의 성공/실패를 추적합니다.
  실패 시 파이프라인을 중단하고 에러 정보를 출력합니다.

  [사용법]
    python run_pipeline.py                # Production (Vina 필수)
    python run_pipeline.py --test-mode    # Test Mode (더미 도킹)
    python run_pipeline.py --skip-boardroom  # Module 03 건너뛰기
    python run_pipeline.py --from 02      # Module 02부터 재시작
    python run_pipeline.py --dashboard    # 파이프라인 완료 후 대시보드 자동 실행

  [환경 변수]
    NUM_SEQUENCES=100000  python run_pipeline.py --test-mode  # 10만 개로 스케일업
    OLLAMA_MODEL=llama3:8b python run_pipeline.py --test-mode  # Ollama 모델 변경
==========================================================================
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

# Windows UTF-8 호환
if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

PROJECT_DIR = Path(__file__).resolve().parent


# ─── 색상 출력 ──────────────────────────────────────────────────────────
def _c(text, color):
    colors = {"green": "32", "red": "31", "yellow": "33", "cyan": "36", "bold": "1"}
    if sys.platform == "win32":
        os.system("")
    code = colors.get(color, "0")
    return f"\033[{code}m{text}\033[0m"


# ─── 모듈 정의 ──────────────────────────────────────────────────────────

MODULES = [
    {
        "id": "00",
        "name": "Module 00: Peptide Generation",
        "script": "00_peptide_gen.py",
        "output": "output/SSTR2_Raw_Library.csv",
        "description": "SSTR2 표적 펩타이드 서열 대규모 생성",
        "extra_args": [],
    },
    {
        "id": "01",
        "name": "Module 01: PK Stability Filter",
        "script": "01_fast_pk_filter.py",
        "output": "output/SSTR2_Filtered_Top100.csv",
        "description": "5대 안정성 지표 기반 Composite Score 필터링",
        "extra_args": [],
    },
    {
        "id": "02",
        "name": "Module 02: 3D Docking + Surrogate Model",
        "script": "02_train_surrogate.py",
        "output": "output/SSTR2_Top100_Boardroom_Ready.csv",
        "description": "1D→3D 변환 + 도킹 + XGBoost 대리 모델",
        "extra_args": [],  # --test-mode는 동적으로 추가
    },
    {
        "id": "03",
        "name": "Module 03: AI Boardroom",
        "script": "03_ai_boardroom.py",
        "output": "output/Final_Winning_Peptides_Report.md",
        "description": "5인 AI 드림팀 원탁회의 평가",
        "extra_args": [],
    },
    {
        "id": "04",
        "name": "Module 04: Validation & Analysis",
        "script": "04_validate_and_analyze.py",
        "output": "output/Validation_Report.md",
        "description": "학술지용 검증, 시각화, K_d 도출",
        "extra_args": [],
    },
]


def run_module(module, test_mode=False):
    """단일 모듈을 실행하고 결과를 반환한다."""
    script_path = PROJECT_DIR / module["script"]
    if not script_path.exists():
        print(f"  {_c('[ERROR]', 'red')} 스크립트 파일 없음: {module['script']}")
        return False

    cmd = [sys.executable, str(script_path)] + module["extra_args"]

    # Module 02에 --test-mode 추가
    if test_mode and module["id"] == "02":
        cmd.append("--test-mode")

    print(f"\n{'='*60}")
    print(f"  {_c(module['name'], 'cyan')}")
    print(f"  {module['description']}")
    print(f"  실행: {' '.join(cmd)}")
    print(f"{'='*60}\n")

    start = time.time()
    try:
        result = subprocess.run(
            cmd,
            cwd=str(PROJECT_DIR),
            env={**os.environ, "PYTHONIOENCODING": "utf-8"},
            timeout=7200,  # 최대 2시간
        )
        elapsed = time.time() - start

        if result.returncode == 0:
            output_path = PROJECT_DIR / module["output"]
            if output_path.exists():
                print(f"\n  {_c('[SUCCESS]', 'green')} {module['name']} "
                      f"완료 ({elapsed:.1f}초)")
                print(f"  출력: {module['output']}")
                return True
            else:
                print(f"\n  {_c('[WARN]', 'yellow')} 프로세스 종료됨, "
                      f"출력 파일 미생성: {module['output']}")
                return True  # 프로세스는 성공했으므로 계속 진행
        else:
            print(f"\n  {_c('[FAIL]', 'red')} {module['name']} "
                  f"실패 (exit code: {result.returncode}, {elapsed:.1f}초)")
            return False

    except subprocess.TimeoutExpired:
        print(f"\n  {_c('[TIMEOUT]', 'red')} {module['name']} — 2시간 초과")
        return False
    except KeyboardInterrupt:
        print(f"\n  {_c('[INTERRUPTED]', 'yellow')} 사용자가 중단했습니다.")
        return False
    except Exception as e:
        print(f"\n  {_c('[ERROR]', 'red')} {module['name']}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="SSTR2 AI Scientist 파이프라인 순차 실행기"
    )
    parser.add_argument("--test-mode", action="store_true",
                        help="Module 02를 더미 도킹 모드로 실행 (Vina 불필요)")
    parser.add_argument("--skip-boardroom", action="store_true",
                        help="Module 03 (AI 원탁회의) 건너뛰기")
    parser.add_argument("--from", dest="start_from", default="00",
                        choices=["00", "01", "02", "03", "04"],
                        help="지정 모듈부터 실행 (이전 출력 파일 필요)")
    parser.add_argument("--dashboard", action="store_true",
                        help="파이프라인 완료 후 대시보드 자동 실행")
    args = parser.parse_args()

    print("\n" + "=" * 60)
    print("  SSTR2 AI Scientist — Pipeline Runner")
    print("=" * 60)
    print(f"  모드        : {'TEST (더미 도킹)' if args.test_mode else 'PRODUCTION'}")
    print(f"  시작 모듈   : Module {args.start_from}")
    print(f"  원탁회의    : {'건너뛰기' if args.skip_boardroom else '실행'}")
    print(f"  대시보드    : {'자동 실행' if args.dashboard else '수동'}")
    print(f"  프로젝트    : {PROJECT_DIR}")
    print("=" * 60)

    if args.test_mode:
        print(f"\n  {_c('WARNING: TEST MODE — 더미 데이터로 실행됩니다.', 'yellow')}")
        print(f"  {_c('학술 출판에 사용할 수 없습니다.', 'yellow')}")

    # 모듈 필터링
    pipeline_start = time.time()
    modules_to_run = [m for m in MODULES if m["id"] >= args.start_from]

    if args.skip_boardroom:
        modules_to_run = [m for m in modules_to_run if m["id"] != "03"]

    success_count = 0
    fail_count = 0

    for module in modules_to_run:
        success = run_module(module, test_mode=args.test_mode)
        if success:
            success_count += 1
        else:
            fail_count += 1
            print(f"\n  {_c('파이프라인 중단:', 'red')} {module['name']}에서 오류 발생")
            print(f"  이전 단계까지의 출력은 output/ 폴더에 보존됩니다.")
            break

    total_elapsed = time.time() - pipeline_start

    # 최종 보고
    print(f"\n{'='*60}")
    print(f"  파이프라인 실행 결과")
    print(f"{'='*60}")
    print(f"  성공: {success_count} / 실패: {fail_count}")
    print(f"  총 소요 시간: {total_elapsed:.1f}초 ({total_elapsed/60:.1f}분)")

    if fail_count == 0:
        print(f"\n  {_c('모든 모듈이 성공적으로 완료되었습니다!', 'green')}")

        # 출력 파일 목록
        print(f"\n  [생성된 파일]")
        output_dir = PROJECT_DIR / "output"
        if output_dir.exists():
            for f in sorted(output_dir.glob("*")):
                if f.is_file():
                    size = f.stat().st_size
                    size_str = (f"{size/1024:.1f}KB" if size < 1024*1024
                                else f"{size/1024/1024:.1f}MB")
                    print(f"    {f.name:<45s} {size_str:>10s}")
            fig_dir = output_dir / "figures"
            if fig_dir.exists():
                n_figs = len(list(fig_dir.glob("*")))
                print(f"    figures/ ({n_figs} files)")
    else:
        print(f"\n  {_c('파이프라인이 중단되었습니다.', 'red')}")

    # 대시보드 자동 실행
    if args.dashboard and fail_count == 0:
        print(f"\n  대시보드를 시작합니다...")
        print(f"  브라우저에서 http://localhost:8501 을 열어주세요.")
        print(f"  종료: Ctrl+C\n")
        subprocess.run(
            [sys.executable, "-m", "streamlit", "run",
             str(PROJECT_DIR / "dashboard.py")],
            cwd=str(PROJECT_DIR),
        )

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
