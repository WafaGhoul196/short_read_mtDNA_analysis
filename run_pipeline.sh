#!/bin/bash
# =============================================================================
#  Mitochondrial DNA Analysis Pipeline
#  Full pipeline: BAM extraction → FASTQ → Alignment & Pileup →
#                 Variant Filtering → Haplogroup Classification →
#                 Haplogroup Variant Removal → Annotation
# =============================================================================

set -uo pipefail

# ── Colours ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

# ── Helpers ───────────────────────────────────────────────────────────────────
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')]${NC} $*" | tee -a "$LOG_FILE"; }
ok()   { echo -e "${GREEN}✔${NC}  $*" | tee -a "$LOG_FILE"; }
warn() { echo -e "${YELLOW}⚠${NC}  $*" | tee -a "$LOG_FILE"; }
die()  { echo -e "${RED}✘  $*${NC}" | tee -a "$LOG_FILE"; exit 1; }

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<EOF
${BOLD}Usage:${NC}
  $0 -c <samples.csv> [OPTIONS]

${BOLD}Required:${NC}
  -c  Path to sample CSV file (must contain a 'Sequencing_number' column at position 6)

${BOLD}Options:${NC}
  -b  Base directory of DRAGEN WGS results  [default: /mnt/DRAGEN_pipeline_results/wgs]
  -o  Output directory                      [default: ./mtdna_results]
  -r  Reference genome FASTA (rCRS)         [default: ./reference/rCRS.fasta]
  -H  Path to Haplogrep3 binary             [default: haplogrep3]
  -t  Threads                               [default: 8]
  -s  Steps to run (comma-separated):
        all | extract | fastq | pileup | filter | haplogrep | haplo_filter | annotate
      [default: all]
  -h  Show this help

${BOLD}Example:${NC}
  $0 -c data_filtered.csv -o results/ -t 16

${BOLD}Dependencies:${NC}
  samtools, bwa, python3 (pandas), R (dplyr, tidyr, ggplot2), haplogrep3
EOF
  exit 0
}

# ── Defaults ──────────────────────────────────────────────────────────────────
BASE_DIR="/mnt/DRAGEN_pipeline_results/wgs"
OUTPUT_DIR="./mtdna_results"
REF_GENOME="./reference/rCRS.fasta"
HAPLOGREP_BIN="haplogrep3"
THREADS=8
STEPS="all"
CSV_FILE=""

# ── Parse arguments ───────────────────────────────────────────────────────────
while getopts "c:b:o:r:H:t:s:h" opt; do
  case $opt in
    c) CSV_FILE="$OPTARG" ;;
    b) BASE_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    r) REF_GENOME="$OPTARG" ;;
    H) HAPLOGREP_BIN="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    s) STEPS="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "$CSV_FILE" ]] && die "Sample CSV (-c) is required. Run '$0 -h' for help."
[[ ! -f "$CSV_FILE" ]] && die "CSV file not found: $CSV_FILE"

# ── Setup directories & log ───────────────────────────────────────────────────
FASTQ_DIR="$OUTPUT_DIR/fastq"
RESULTS_DIR="$OUTPUT_DIR/results"
HAPLO_DIR="$OUTPUT_DIR/haplogrep"
HAPLO_FILTER_DIR="$OUTPUT_DIR/haplogroup_filtered"
ANNOT_DIR="$OUTPUT_DIR/annotation"
LOG_FILE="$OUTPUT_DIR/pipeline.log"

mkdir -p "$OUTPUT_DIR" "$FASTQ_DIR" "$RESULTS_DIR" "$HAPLO_DIR" "$HAPLO_FILTER_DIR" "$ANNOT_DIR"

# Resolve script directory (for helper scripts)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "" | tee -a "$LOG_FILE"
echo -e "${BOLD}╔══════════════════════════════════════════════════╗${NC}" | tee -a "$LOG_FILE"
echo -e "${BOLD}║     Mitochondrial DNA Analysis Pipeline          ║${NC}" | tee -a "$LOG_FILE"
echo -e "${BOLD}╚══════════════════════════════════════════════════╝${NC}" | tee -a "$LOG_FILE"
log "Pipeline started"
log "CSV        : $CSV_FILE"
log "Output     : $OUTPUT_DIR"
log "Reference  : $REF_GENOME"
log "Threads    : $THREADS"
log "Steps      : $STEPS"

# ── Dependency check ──────────────────────────────────────────────────────────
check_deps() {
  log "Checking dependencies..."
  local missing=()
  for cmd in samtools bwa python3 Rscript; do
    command -v "$cmd" &>/dev/null || missing+=("$cmd")
  done
  if [[ "$STEPS" == *"haplogrep"* || "$STEPS" == "all" ]]; then
    command -v "$HAPLOGREP_BIN" &>/dev/null || missing+=("$HAPLOGREP_BIN")
  fi
  if [[ ${#missing[@]} -gt 0 ]]; then
    die "Missing dependencies: ${missing[*]}\nSee README.md for installation instructions."
  fi
  ok "All dependencies found."
}

run_step() { [[ "$STEPS" == "all" || "$STEPS" == *"$1"* ]]; }

# ── Read sample IDs ───────────────────────────────────────────────────────────
mapfile -t SAMPLE_IDS < <(tail -n +2 "$CSV_FILE" | cut -d',' -f6 | xargs -I{} echo {})

# =============================================================================
# STEP 1 – Extract mitochondrial BAM
# =============================================================================
step_extract() {
  log "━━━ STEP 1: Extract mitochondrial BAM ━━━"
  local success=0 errors=0

  for raw_id in "${SAMPLE_IDS[@]}"; do
    local seq_id
    seq_id=$(echo "$raw_id" | xargs)
    local bam_file="$BASE_DIR/$seq_id/dragen/${seq_id}.bam"
    local output_bam="$OUTPUT_DIR/${seq_id}_chrM.bam"

    if [[ -f "$output_bam" && -f "${output_bam}.bai" ]]; then
      warn "$seq_id: chrM BAM already exists, skipping."
      success=$((success + 1)); continue
    fi

    if [[ ! -d "$BASE_DIR/$seq_id/dragen" ]]; then
      warn "$seq_id: dragen folder not found."; errors=$((errors + 1)); continue
    fi
    if [[ ! -f "$bam_file" ]]; then
      warn "$seq_id: BAM file not found at $bam_file"; errors=$((errors + 1)); continue
    fi

    local chr
    if samtools view -H "$bam_file" | grep -q -E "@SQ\s+SN:chrM"; then chr="chrM"
    elif samtools view -H "$bam_file" | grep -q -E "@SQ\s+SN:MT"; then chr="MT"
    else warn "$seq_id: mitochondrial chromosome not found."; errors=$((errors + 1)); continue; fi

    log "$seq_id: extracting $chr..."
    samtools view -h "$bam_file" "$chr" | samtools view -b -o "$output_bam" 2>>"$LOG_FILE"
    samtools index "$output_bam" 2>>"$LOG_FILE"
    ok "$seq_id: extracted ($chr)"; success=$((success + 1))
  done

  log "Step 1 complete — success: $success, errors: $errors"
  if [[ $errors -gt 0 ]]; then warn "$errors samples had errors. Check $LOG_FILE."; fi
}

# =============================================================================
# STEP 2 – BAM to FASTQ
# =============================================================================
step_fastq() {
  log "━━━ STEP 2: BAM → FASTQ ━━━"

  for bam_file in "$OUTPUT_DIR"/*_chrM.bam; do
    [[ -e "$bam_file" ]] || { warn "No chrM BAM files found in $OUTPUT_DIR"; return; }
    local base_name
    base_name=$(basename "$bam_file" .bam)
    local r1="$FASTQ_DIR/${base_name}_R1.fq.gz"

    if [[ -f "$r1" ]]; then
      warn "$base_name: FASTQ already exists, skipping."; continue
    fi

    log "$base_name: sorting by name..."
    local sorted="$FASTQ_DIR/${base_name}.name_sorted.bam"
    samtools sort -n -@ "$THREADS" -o "$sorted" "$bam_file"

    log "$base_name: converting to FASTQ..."
    samtools fastq -@ "$THREADS" \
      -1 "$FASTQ_DIR/${base_name}_R1.fq.gz" \
      -2 "$FASTQ_DIR/${base_name}_R2.fq.gz" \
      -s "$FASTQ_DIR/${base_name}_unpaired.fq.gz" \
      "$sorted"

    rm -f "$sorted"
    ok "$base_name: FASTQ generated."
  done
}

# =============================================================================
# STEP 3 – Align & generate pileup
# =============================================================================
step_pileup() {
  log "━━━ STEP 3: Alignment & Pileup ━━━"
  [[ ! -f "$REF_GENOME" ]] && die "Reference genome not found: $REF_GENOME"

  # Index reference if needed
  if [[ ! -f "${REF_GENOME}.bwt" ]]; then
    log "Indexing reference genome..."
    bwa index "$REF_GENOME" 2>>"$LOG_FILE"
    samtools faidx "$REF_GENOME" 2>>"$LOG_FILE"
    ok "Reference indexed."
  fi

  for READ1 in "$FASTQ_DIR"/*_R1.fq.gz; do
    [[ -e "$READ1" ]] || { warn "No FASTQ R1 files found."; return; }
    local SAMPLE
    SAMPLE=$(basename "$READ1" _R1.fq.gz)
    local READ2="$FASTQ_DIR/${SAMPLE}_R2.fq.gz"
    local MITO_BAM="$RESULTS_DIR/mitochondrial_${SAMPLE}.bam"
    local PILEUP_FILE="$RESULTS_DIR/mitochondrial_pileup_${SAMPLE}.txt"
    # Output directly with _heteroplasmy suffix — no separate recalculation step needed
    local OUTPUT_CSV="$RESULTS_DIR/output_pileup_analysis_${SAMPLE}_heteroplasmy.csv"

    if [[ -f "$OUTPUT_CSV" ]]; then
      warn "$SAMPLE: pileup CSV already exists, skipping."; continue
    fi

    log "$SAMPLE: aligning with BWA MEM..."
    local SORTED_BAM="$RESULTS_DIR/tmp_${SAMPLE}_sorted.bam"

    bwa mem -t "$THREADS" "$REF_GENOME" "$READ1" "$READ2" 2>>"$LOG_FILE" \
      | samtools view -@ "$THREADS" -bS -h - \
      | samtools sort -@ "$THREADS" -o "$SORTED_BAM" -
    samtools index "$SORTED_BAM"

    log "$SAMPLE: extracting mtDNA reads..."
    local MITO_CONTIG
    if samtools view -H "$SORTED_BAM" | grep -q "SN:NC_012920"; then
      MITO_CONTIG="NC_012920.1"
    elif samtools view -H "$SORTED_BAM" | grep -q "SN:chrM"; then
      MITO_CONTIG="chrM"
    elif samtools view -H "$SORTED_BAM" | grep -q "SN:MT"; then
      MITO_CONTIG="MT"
    else
      warn "$SAMPLE: mitochondrial contig not detected; using full BAM."; MITO_CONTIG=""
    fi

    if [[ -n "$MITO_CONTIG" ]]; then
      samtools view -@ "$THREADS" -b -h -F 4 "$SORTED_BAM" "$MITO_CONTIG" > "$MITO_BAM"
    else
      samtools view -@ "$THREADS" -b -h -F 4 "$SORTED_BAM" > "$MITO_BAM"
    fi
    samtools index "$MITO_BAM"
    rm -f "$SORTED_BAM" "${SORTED_BAM}.bai"

    log "$SAMPLE: generating pileup..."
    samtools mpileup -f "$REF_GENOME" -o "$PILEUP_FILE" \
      -A -d 15000 -Q 30 -B "$MITO_BAM" 2>>"$LOG_FILE"

    log "$SAMPLE: parsing pileup to CSV (heteroplasmy calculated inline)..."
    python3 "$SCRIPT_DIR/scripts/pileup_analysis.py" "$PILEUP_FILE" "$OUTPUT_CSV"

    ok "$SAMPLE: pileup complete → $OUTPUT_CSV"
  done
}

# =============================================================================
# STEP 4 – Filter variants (R)
# =============================================================================
step_filter() {
  log "━━━ STEP 4: Filter Variants (R) ━━━"
  Rscript "$SCRIPT_DIR/scripts/pileup_filter.R" "$RESULTS_DIR" 2>>"$LOG_FILE"
  ok "Variant filtering complete."
}

# =============================================================================
# STEP 5 – Haplogroup classification (Haplogrep3)
# =============================================================================
step_haplogrep() {
  log "━━━ STEP 5: Haplogroup Classification ━━━"
  local OUTPUT="$HAPLO_DIR/merged_haplogrep_results.txt"
  echo -e "SampleID\tHaplogroup\tRank\tQuality\tRange\tNot_Found_Polys\tFound_Polys\tRemaining_Polys\tAAC_In_Remainings\tInput_Sample" > "$OUTPUT"

  for csv in "$RESULTS_DIR"/*_data_filtred_heteroplasmy.csv; do
    [[ -e "$csv" ]] || { warn "No filtered CSV files found in $RESULTS_DIR"; return; }
    local sample
    sample=$(basename "$csv" _data_filtred_heteroplasmy.csv)
    log "$sample: converting CSV → VCF..."

    local raw_vcf="$HAPLO_DIR/${sample}_raw.vcf"
    local fmt_vcf="$HAPLO_DIR/${sample}_formatted.vcf"

    python3 "$SCRIPT_DIR/scripts/csv_to_vcf.py" -i "$csv" -o "$raw_vcf"
    [[ ! -f "$raw_vcf" ]] && { warn "$sample: VCF conversion failed."; continue; }

    awk -v sample="$sample" '
    BEGIN {
      FS=OFS="\t";
      print "##fileformat=VCFv4.2";
      print "##FILTER=<ID=PASS,Description=\"Variants passed filters\">";
      print "##FORMAT=<ID=AF,Number=1,Type=String,Description=\"Allele Frequency\">";
      print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">";
      print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
      print "##contig=<ID=NC_012920.1,length=16569>";
      print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" sample
    }
    !/^#/ {
      split($8, info, ";");
      for (i in info) {
        split(info[i], arr, "=");
        if (arr[1] == "AF") af = arr[2];
        if (arr[1] == "DP") dp = arr[2];
      }
      gt = (af > 0.9) ? "1/1" : "0/1";
      print "NC_012920.1", $2, ".", $4, $5, ".", "PASS", ".", "GT:AF:DP", gt ":" af ":" dp
    }' "$raw_vcf" > "$fmt_vcf"
    rm -f "$raw_vcf"

    log "$sample: running Haplogrep3..."
    local haplo_out="$HAPLO_DIR/${sample}_haplogrep.txt"
    "$HAPLOGREP_BIN" classify \
      --tree phylotree-rcrs@17.2 \
      --input "$fmt_vcf" \
      --output "$haplo_out" \
      --extend-report 2>>"$LOG_FILE"

    if [[ -f "$haplo_out" ]]; then
      tail -n +2 "$haplo_out" >> "$OUTPUT"
      rm -f "$haplo_out"
      ok "$sample: haplogroup classified."
    else
      warn "$sample: Haplogrep failed."
    fi
  done
  ok "Haplogrep complete → $OUTPUT"
}

# =============================================================================
# STEP 6 – Remove haplogroup-defining variants
# =============================================================================
step_haplo_filter() {
  log "━━━ STEP 6: Filter Haplogroup-Defining Variants ━━━"
  local haplo_results="$HAPLO_DIR/merged_haplogrep_results.txt"
  if [[ ! -f "$haplo_results" ]]; then
    die "Haplogrep results not found: $haplo_results\nRun the haplogrep step first."
  fi
  python3 "$SCRIPT_DIR/scripts/filter_haplo_variants.py" \
    --haplogrep  "$haplo_results" \
    --input_dir  "$RESULTS_DIR" \
    --output_dir "$HAPLO_FILTER_DIR" \
    2>>"$LOG_FILE"
  ok "Haplogroup variant filtering complete → $HAPLO_FILTER_DIR"
}

# =============================================================================
# STEP 7 – Annotation
# =============================================================================
step_annotate() {
  log "━━━ STEP 7: Annotation ━━━"
  python3 "$SCRIPT_DIR/scripts/apply_annotation.py" \
    --input_dir  "$HAPLO_FILTER_DIR" \
    --annot_dir  "$SCRIPT_DIR/annotation_databases/curated" \  # either curated or raw folder depending on user choice
    --output     "$ANNOT_DIR/merged_annotated.csv" \
    2>>"$LOG_FILE"
  ok "Annotation complete → $ANNOT_DIR/merged_annotated.csv"
}

# =============================================================================
# MAIN
# =============================================================================
check_deps

run_step "extract"      && step_extract      || true
run_step "fastq"        && step_fastq        || true
run_step "pileup"       && step_pileup       || true
run_step "filter"       && step_filter       || true
run_step "haplogrep"    && step_haplogrep    || true
run_step "haplo_filter" && step_haplo_filter || true
run_step "annotate"     && step_annotate     || true

echo "" | tee -a "$LOG_FILE"
echo -e "${GREEN}${BOLD}Pipeline completed successfully at $(date)${NC}" | tee -a "$LOG_FILE"
echo -e "Results are in: ${BOLD}$OUTPUT_DIR${NC}" | tee -a "$LOG_FILE"
