#!/usr/bin/env bash
# =============================================================================
#  GermCall test on a tiny SIMULATED data set (a few minutes). Needs all tools of
#  environment.yml (incl. GATK). Output goes to ./germcall_test/ (safe to delete).
#      bash examples/run_test.sh
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$(pwd)/germcall_test"; mkdir -p "${OUT}"
python3 "${ROOT}/examples/make_test_data.py" "${OUT}/data"
cat > "${OUT}/test.yaml" <<YAML
project: {name: test, species: "Simulated plant", prefix: test, workdir: ${OUT}/work, threads: 4}
calling:
  reference: ${OUT}/data/ref.fa
  samples_file: ${OUT}/data/samples.tsv
  paired: true
  parallel_samples: 2
  threads_per_sample: 2
  parallel_intervals: 2
  java_mem_gb: 4
  n_chromosomes: 2
YAML
bash "${ROOT}/germcall.sh" --config "${OUT}/test.yaml" --all
echo "Test finished. Report: ${OUT}/work/qc_reports/multiqc_report.html"
