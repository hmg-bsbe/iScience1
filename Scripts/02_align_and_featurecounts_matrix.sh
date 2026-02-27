#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Multi-dataset paired-end RNA-seq pipeline:
#   (Optional) strandedness inference (Salmon) 
#   HISAT2 indexing (one-time) 
#   HISAT2 alignment (per sample) 
#   featureCounts (single combined matrix) 
#   counts-only matrix export (GeneID + raw counts)
#
# Usage:
#   bash align_and_featurecounts_matrix.sh samples.tsv
#
# samples.tsv (TAB-separated, with header):
#   sample_id    read1                           read2
#   S1           /path/S1_R1.fastq.gz             /path/S1_R2.fastq.gz
#   ...
# ------------------------------------------------------------

THREADS=8
OUTDIR="./align_counts_out"

# --------- REQUIRED: edit these reference paths ----------
GENOME_FA="/path/to/genome.fa"            # we used Mus musculus GRCm39 genome FASTA
GTF="/path/to/annotation.gtf"             # matching the same genome build
HISAT2_INDEX_PREFIX="/path/to/hisat2_index/genome_prefix"

# Optional strandedness inference (Salmon):
TX_FASTA="/path/to/transcripts.cdna.fa"   # transcriptome FASTA (cDNA/transcripts)
SALMON_INDEX="/path/to/salmon_tx_index"   # Salmon index folder
# --------------------------------------------------------

# Strandedness control:
#   "salmon" -> infer per sample (stores a summary table; alignment uses inferred settings)
#   "manual" -> use fixed settings below across all samples
STRANDNESS_MODE="salmon"   # "salmon" or "manual"

# Manual settings (used only if STRANDNESS_MODE="manual")
HISAT2_STRAND_MANUAL="RF"  # FR or RF (paired-end)
FC_STRAND_MANUAL=2         # 0 unstranded, 1 forward, 2 reverse

SAMPLES_TSV="${1:?Provide samples.tsv}"

mkdir -p "${OUTDIR}"/{strandness,hisat2,bam,counts,logs,tmp}

# ---- Tool checks ----
for tool in hisat2 hisat2-build samtools featureCounts; do
  command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not found in PATH"; exit 1; }
done

if [[ "${STRANDNESS_MODE}" == "salmon" ]]; then
  command -v salmon >/dev/null 2>&1 || { echo "ERROR: salmon not found in PATH (needed for STRANDNESS_MODE=salmon)"; exit 1; }
  [[ -f "${TX_FASTA}" ]] || { echo "ERROR: TX_FASTA not found: ${TX_FASTA}"; exit 1; }
fi

[[ -f "${GENOME_FA}" ]] || { echo "ERROR: GENOME_FA not found: ${GENOME_FA}"; exit 1; }
[[ -f "${GTF}" ]]      || { echo "ERROR: GTF not found: ${GTF}"; exit 1; }
[[ -f "${SAMPLES_TSV}" ]] || { echo "ERROR: samples.tsv not found: ${SAMPLES_TSV}"; exit 1; }

# ---- Salmon index (one-time, only if using salmon mode) ----
if [[ "${STRANDNESS_MODE}" == "salmon" ]]; then
  echo "[0/5] Preparing Salmon index (one-time) if missing..."
  if [[ ! -d "${SALMON_INDEX}" ]]; then
    mkdir -p "$(dirname "${SALMON_INDEX}")"
    salmon index -t "${TX_FASTA}" -i "${SALMON_INDEX}" -k 31
  else
    echo "  Salmon index exists: ${SALMON_INDEX}"
  fi
fi

# ---- HISAT2 index (one-time) ----
echo "[1/5] Preparing HISAT2 index (one-time) if missing..."
if ls "${HISAT2_INDEX_PREFIX}".*.ht2 >/dev/null 2>&1; then
  echo "  HISAT2 index exists: ${HISAT2_INDEX_PREFIX}.*.ht2"
else
  mkdir -p "$(dirname "${HISAT2_INDEX_PREFIX}")"
  hisat2-build -p "${THREADS}" "${GENOME_FA}" "${HISAT2_INDEX_PREFIX}"
fi

# ---- Helpers ----
extract_salmon_libtype() {
  local log_file="$1"
  [[ -f "${log_file}" ]] || { echo "NA"; return; }
  # line: "Automatically detected most likely library type as ISR"
  local lib
  lib="$(grep -i "most likely library type" "${log_file}" | tail -n 1 | awk '{print $NF}' || true)"
  [[ -n "${lib}" ]] && echo "${lib}" || echo "NA"
}

map_salmon_to_settings() {
  # Maps Salmon libtype to (HISAT2_STRAND, FC_STRAND)
  # ISR / OSR -> reverse stranded paired-end => HISAT2 RF, featureCounts -s 2
  # ISF / OSF -> forward stranded paired-end => HISAT2 FR, featureCounts -s 1
  # otherwise -> unknown/unstranded => HISAT2 no strand flag, featureCounts -s 0
  local libtype="$1"
  case "${libtype}" in
    ISR|OSR) echo -e "RF\t2" ;;
    ISF|OSF) echo -e "FR\t1" ;;
    *)       echo -e "NA\t0" ;;
  esac
}

# ---- Summary file ----
echo -e "sample_id\tsalmon_libtype\thisat2_strand\tfeatureCounts_s" > "${OUTDIR}/logs/strandness_settings.tsv"
: > "${OUTDIR}/tmp/bam_list.txt"

# ---- Per-sample alignment ----
echo "[2/5] Aligning samples..."
tail -n +2 "${SAMPLES_TSV}" | while IFS=$'\t' read -r SAMPLE_ID READ1 READ2; do
  [[ -n "${SAMPLE_ID}" ]] || continue
  [[ -f "${READ1}" ]] || { echo "ERROR: Missing READ1 for ${SAMPLE_ID}: ${READ1}"; exit 1; }
  [[ -f "${READ2}" ]] || { echo "ERROR: Missing READ2 for ${SAMPLE_ID}: ${READ2}"; exit 1; }

  echo "  → ${SAMPLE_ID}"

  SALMON_LIBTYPE="NA"
  HISAT2_STRAND="${HISAT2_STRAND_MANUAL}"
  FC_STRAND="${FC_STRAND_MANUAL}"

  if [[ "${STRANDNESS_MODE}" == "salmon" ]]; then
    echo "    [Salmon] inferring strandedness..."
    SALMON_OUT="${OUTDIR}/strandness/${SAMPLE_ID}"
    salmon quant \
      -i "${SALMON_INDEX}" \
      -l A \
      -1 "${READ1}" \
      -2 "${READ2}" \
      -p "${THREADS}" \
      --validateMappings \
      -o "${SALMON_OUT}" \
      > "${OUTDIR}/logs/${SAMPLE_ID}.salmon.stdout.log" \
      2> "${OUTDIR}/logs/${SAMPLE_ID}.salmon.stderr.log"

    SALMON_LOG="${SALMON_OUT}/logs/salmon_quant.log"
    SALMON_LIBTYPE="$(extract_salmon_libtype "${SALMON_LOG}")"

    read -r HISAT2_STRAND FC_STRAND < <(map_salmon_to_settings "${SALMON_LIBTYPE}")
  fi

  echo -e "${SAMPLE_ID}\t${SALMON_LIBTYPE}\t${HISAT2_STRAND}\t${FC_STRAND}" >> "${OUTDIR}/logs/strandness_settings.tsv"

  echo "    [HISAT2] aligning and sorting..."
  BAM_OUT="${OUTDIR}/bam/${SAMPLE_ID}.sorted.bam"
  HISAT2_LOG="${OUTDIR}/hisat2/${SAMPLE_ID}.hisat2.log"

  if [[ "${HISAT2_STRAND}" == "FR" || "${HISAT2_STRAND}" == "RF" ]]; then
    hisat2 -p "${THREADS}" \
      --rna-strandness "${HISAT2_STRAND}" \
      -x "${HISAT2_INDEX_PREFIX}" \
      -1 "${READ1}" -2 "${READ2}" \
      2> "${HISAT2_LOG}" \
    | samtools view -@ "${THREADS}" -bS - \
    | samtools sort -@ "${THREADS}" -o "${BAM_OUT}" -
  else
    # Unknown/unstranded: run without --rna-strandness
    hisat2 -p "${THREADS}" \
      -x "${HISAT2_INDEX_PREFIX}" \
      -1 "${READ1}" -2 "${READ2}" \
      2> "${HISAT2_LOG}" \
    | samtools view -@ "${THREADS}" -bS - \
    | samtools sort -@ "${THREADS}" -o "${BAM_OUT}" -
  fi

  samtools index "${BAM_OUT}"

  # Track BAMs for a single combined featureCounts run
  echo "${BAM_OUT}" >> "${OUTDIR}/tmp/bam_list.txt"
done

# ---- Decide one featureCounts strandedness setting for the combined matrix ----
# If STRANDNESS_MODE=salmon and samples disagree, we default to unstranded (0) to avoid forcing a wrong setting.
echo "[3/5] Determining featureCounts strandedness for the combined run..."
FC_STRAND_FINAL="${FC_STRAND_MANUAL}"

if [[ "${STRANDNESS_MODE}" == "salmon" ]]; then
  # Collect unique featureCounts_s values observed (excluding header)
  UNIQUE_S="$(tail -n +2 "${OUTDIR}/logs/strandness_settings.tsv" | cut -f4 | sort -u | tr '\n' ' ')"
  # If exactly one unique value (0 or 1 or 2), use it; else use 0.
  COUNT_UNIQUE="$(tail -n +2 "${OUTDIR}/logs/strandness_settings.tsv" | cut -f4 | sort -u | wc -l | awk '{print $1}')"
  if [[ "${COUNT_UNIQUE}" -eq 1 ]]; then
    FC_STRAND_FINAL="$(tail -n +2 "${OUTDIR}/logs/strandness_settings.tsv" | cut -f4 | head -n 1)"
  else
    FC_STRAND_FINAL=0
    echo "  NOTE: Mixed/uncertain strandedness across samples; using featureCounts -s 0 (unstranded) for the combined matrix."
    echo "  You can split datasets by model/batch and rerun counts if needed."
  fi
fi

echo "  featureCounts -s ${FC_STRAND_FINAL}"

# ---- Run featureCounts ONCE on all BAMs (combined matrix) ----
echo "[4/5] Running featureCounts (combined matrix)..."
COUNTS_FULL="${OUTDIR}/counts/featureCounts.full.txt"
COUNTS_LOG="${OUTDIR}/logs/featureCounts.combined.log"

# Read BAM list into an array
mapfile -t BAMS < "${OUTDIR}/tmp/bam_list.txt"
if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "ERROR: No BAMs found in ${OUTDIR}/tmp/bam_list.txt"
  exit 1
fi

featureCounts -T "${THREADS}" \
  -s "${FC_STRAND_FINAL}" \
  -a "${GTF}" \
  -o "${COUNTS_FULL}" \
  "${BAMS[@]}" \
  > "${COUNTS_LOG}" 2>&1

# ---- Export counts-only matrix (GeneID + raw counts) ----
echo "[5/5] Exporting GeneID + raw counts matrix..."
COUNTS_CLEAN="${OUTDIR}/counts/featureCounts.gene_counts.txt"
cut -f1,7- "${COUNTS_FULL}" > "${COUNTS_CLEAN}"

echo "[DONE]"
echo "Full featureCounts table: ${COUNTS_FULL}"
echo "Counts-only matrix:       ${COUNTS_CLEAN}"
echo "Strandness summary:       ${OUTDIR}/logs/strandness_settings.tsv"
echo "featureCounts log:        ${COUNTS_LOG}"
