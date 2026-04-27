#!/bin/bash
# =============================================================================
#  mtDNA Pipeline - One-time Setup Script
#  Run this once before using the pipeline for the first time.
#  Usage: bash setup.sh
# =============================================================================

set -uo pipefail

# ── Colours ───────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'
ok()   { echo -e "${GREEN}✔${NC}  $*"; }
warn() { echo -e "${YELLOW}⚠${NC}  $*"; }
die()  { echo -e "${RED}✘  $*${NC}"; exit 1; }
log()  { echo -e "${BOLD}──────────────────────────────────────────${NC}"; echo -e "  $*"; }

echo ""
echo -e "${BOLD}╔══════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}║        mtDNA Pipeline - Setup                    ║${NC}"
echo -e "${BOLD}╚══════════════════════════════════════════════════╝${NC}"
echo ""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
# STEP 1 – Conda environment
# =============================================================================
log "STEP 1: Conda environment"

if ! command -v conda &>/dev/null; then
    die "conda not found. Please install Miniconda first:\n  https://docs.conda.io/en/latest/miniconda.html"
fi

if conda env list | grep -q "mtdna_pipeline"; then
    warn "Conda environment 'mtdna_pipeline' already exists — skipping creation."
    warn "To recreate it from scratch: conda env remove -n mtdna_pipeline"
else
    echo "Creating conda environment from environment.yml..."
    # Ensure correct channel priority before creating
    conda config --add channels conda-forge 2>/dev/null || true
    conda config --set channel_priority strict 2>/dev/null || true
    conda env create -f "$SCRIPT_DIR/environment.yml" \
        || die "Conda environment creation failed. Check environment.yml."
    ok "Conda environment 'mtdna_pipeline' created."
fi

# =============================================================================
# STEP 2 – Haplogrep3
# =============================================================================
log "STEP 2: Haplogrep3"

HAPLOGREP_VERSION="3.2.2"
HAPLOGREP_DIR="/opt/haplogrep3"
HAPLOGREP_BIN="/usr/local/bin/haplogrep3"
HAPLOGREP_ZIP="$SCRIPT_DIR/haplogrep3-${HAPLOGREP_VERSION}-linux.zip"

# Check if already installed
if [[ -f "$HAPLOGREP_BIN" ]] && "$HAPLOGREP_BIN" 2>&1 | grep -q "Haplogrep"; then
    ok "Haplogrep3 already installed at $HAPLOGREP_BIN"
else
    # Check Java
    if ! command -v java &>/dev/null; then
        warn "Java not found — installing via conda..."
        conda run -n mtdna_pipeline conda install -y -c conda-forge openjdk=17 \
            || die "Java installation failed. Install Java 11+ manually."
        ok "Java installed."
    else
        JAVA_VER=$(java -version 2>&1 | head -1)
        ok "Java found: $JAVA_VER"
    fi

    # Download Haplogrep3 if not already downloaded
    if [[ ! -f "$HAPLOGREP_ZIP" ]]; then
        echo "Downloading Haplogrep3 v${HAPLOGREP_VERSION}..."
        wget -q --show-progress \
            "https://github.com/genepi/haplogrep3/releases/download/v${HAPLOGREP_VERSION}/haplogrep3-${HAPLOGREP_VERSION}-linux.zip" \
            -O "$HAPLOGREP_ZIP" \
            || die "Download failed. Check your internet connection."
        ok "Downloaded: $HAPLOGREP_ZIP"
    else
        ok "Zip already downloaded: $HAPLOGREP_ZIP"
    fi

    # Install
    echo "Installing Haplogrep3 to $HAPLOGREP_DIR..."
    TMP_DIR=$(mktemp -d)
    unzip -q "$HAPLOGREP_ZIP" -d "$TMP_DIR"

    sudo mkdir -p "$HAPLOGREP_DIR"
    sudo cp "$TMP_DIR/haplogrep3.jar"  "$HAPLOGREP_DIR/"
    sudo cp "$TMP_DIR/haplogrep3.yaml" "$HAPLOGREP_DIR/" 2>/dev/null || true
    rm -rf "$TMP_DIR"

    # Create launcher script
    sudo tee "$HAPLOGREP_DIR/haplogrep3" > /dev/null << 'EOF'
#!/bin/bash
java -jar /opt/haplogrep3/haplogrep3.jar "$@"
EOF
    sudo chmod +x "$HAPLOGREP_DIR/haplogrep3"

    # Symlink to PATH
    sudo ln -sf "$HAPLOGREP_DIR/haplogrep3" "$HAPLOGREP_BIN"

    ok "Haplogrep3 installed."
fi

# Verify
HAPLO_OUT=$("$HAPLOGREP_BIN" 2>&1 | head -2)
ok "Haplogrep3 check: $(echo "$HAPLO_OUT" | head -1)"

# =============================================================================
# STEP 3 – Reference genome
# =============================================================================
log "STEP 3: Reference genome"

REF="$SCRIPT_DIR/reference/sequence.fasta"
mkdir -p "$SCRIPT_DIR/reference"

if [[ ! -f "$REF" ]]; then
    warn "Reference genome not found at: $REF"
    echo ""
    echo "  Please copy your reference genome there before running the pipeline:"
    echo "    cp /path/to/your/reference.fasta $REF"
    echo ""
    echo "  Or use a different path with the -r flag when running run_pipeline.sh."
    echo ""
else
    ok "Reference genome found: $REF"

    # Index if needed
    if [[ ! -f "${REF}.bwt" ]]; then
        echo "Indexing reference genome (BWA + samtools)..."
        conda run -n mtdna_pipeline bash -c "
            bwa index '$REF' && samtools faidx '$REF'
        " || warn "Indexing failed — will be retried automatically during the first pipeline run."
        ok "Reference indexed."
    else
        ok "Reference already indexed."
    fi
fi

# =============================================================================
# Done
# =============================================================================
echo ""
echo -e "${GREEN}${BOLD}╔══════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}${BOLD}║  Setup complete! Ready to run the pipeline.      ║${NC}"
echo -e "${GREEN}${BOLD}╚══════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${BOLD}1) Activate the environment:${NC}"
echo -e "   conda activate mtdna_pipeline"
echo ""

echo -e "${BOLD}2) Run the pipeline:${NC}"
echo ""

echo "bash run_pipeline.sh \\"
echo "  -c <samples.csv> \\"
echo "  -b <DRAGEN_results_dir> \\"
echo "  -o <output_dir> \\"
echo "  -r reference/rCRS.fasta \\"
echo "  -t 16"
echo ""

echo -e "${BOLD}Optional parameters:${NC}"
echo ""
echo "  -H /usr/local/bin/haplogrep3   # Haplogrep path"
echo "  -s <steps>                     # Run specific steps"
echo "  -n false                       # Skip annotation"
echo ""
echo -e "${BOLD}Filtering parameters (Step 4):${NC}"
echo ""
echo "  -m 0.03   # Minimum heteroplasmy"
echo "  -a 0.03   # Minimum allele frequency (with heteroplasmy)"
echo "  -A 0.05   # Strict allele frequency threshold"
echo "  -u 0.5    # Max unmapped/depth ratio"
echo ""

echo -e "${BOLD}Example:${NC}"
echo ""
echo "bash run_pipeline.sh \\"
echo "  -c samples.csv \\"
echo "  -o results/ \\"
echo "  -t 16 \\"
echo "  -m 0.03 -a 0.03 -A 0.05 -u 0.5 \\"
echo "  -n false"
echo ""

echo -e "${CYAN}Tip:${NC} You can run specific steps using -s (e.g. haplogrep,annotate)"
echo ""
