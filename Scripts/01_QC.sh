#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# QC: FastQC + MultiQC
# -----------------------------

# Arguments:
#   1) FASTQ_DIR  : directory containing .fastq.gz files
#   2) OUT_DIR    : output directory for QC results
#   3) THREADS    : threads for FastQC (default: 4)

FASTQ_DIR="${1:?Provide FASTQ input directory}"
OUT_DIR="${2:-qc_out}"
THREADS="${3:-4}"

FASTQC_OUT="${OUT_DIR}/fastqc"
MULTIQC_OUT="${OUT_DIR}/multiqc"

mkdir -p "${FASTQC_OUT}" "${MULTIQC_OUT}"

# Check tools
command -v fastqc >/dev/null 2>&1 || { echo "ERROR: fastqc not found in PATH"; exit 1; }
command -v multiqc >/dev/null 2>&1 || { echo "ERROR: multiqc not found in PATH"; exit 1; }

# Check input
shopt -s nullglob
FASTQS=( "${FASTQ_DIR}"/*.fastq.gz "${FASTQ_DIR}"/*.fq.gz )
if [[ ${#FASTQS[@]} -eq 0 ]]; then
  echo "ERROR: No .fastq.gz/.fq.gz files found in ${FASTQ_DIR}"
  exit 1
fi

echo "[INFO] Running FastQC on ${#FASTQS[@]} files with ${THREADS} threads..."
fastqc -t "${THREADS}" -o "${FASTQC_OUT}" "${FASTQS[@]}"

echo "[INFO] Aggregating QC with MultiQC..."
multiqc "${FASTQC_OUT}" -o "${MULTIQC_OUT}" --force

echo "[DONE] FastQC reports:   ${FASTQC_OUT}"
echo "[DONE] MultiQC report:   ${MULTIQC_OUT}/multiqc_report.html"
