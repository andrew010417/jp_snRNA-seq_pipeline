#!/usr/bin/env bash
# 전체 파이프라인을 순서대로 실행합니다.
#   사용법:  bash scripts/run_all.sh [Rscript 경로]
# 각 단계는 체크포인트를 남기므로, 중간에 멈춰도 다시 실행하면 이어서 갑니다.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
RS="${1:-Rscript}"

export SNRNA_ROOT="$ROOT"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export OMP_NUM_THREADS="$OPENBLAS_NUM_THREADS"

mkdir -p "$ROOT/run/logs"

for step in 01_load_qc 02_doublet 03_filter 04_normalize 05_hvg_scale_pca 06_integrate 07_cluster 08_markers; do
  echo ""
  echo "=============================================================="
  echo "  $step"
  echo "=============================================================="
  "$RS" --vanilla "$HERE/${step}.R" 2>&1 | tee -a "$ROOT/run/logs/${step}.log"
done

echo ""
echo "완료. 결과는 $ROOT/run/ 아래에 있습니다."
