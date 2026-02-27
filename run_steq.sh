#!/usr/bin/env bash
# ==============================================================================
# run_steq.sh — Build (if needed) and run STEQ
#
# STEQ (Species Tree Estimation using Quartet) infers a species tree from
# a collection of gene trees in Newick format.  This wrapper handles
# compilation, argument validation, run-time logging, and timing.
# ==============================================================================
set -euo pipefail

# ── Locate project directories ───────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUMMARIZER_DIR="$SCRIPT_DIR/Summarizer"
BINARIES_DIR="$SCRIPT_DIR/Binaries"
STEQ_EXEC="$SUMMARIZER_DIR/STEQ.out"
FASTME_EXEC="$BINARIES_DIR/fastme-2.1.5-linux64"
BIONJ_EXEC="$BINARIES_DIR/BIONJ"

# ── Defaults ─────────────────────────────────────────────────────────────────
INPUT=""
OUTPUT=""
MATRIX=""
CENTRAL_TENDENCY=0
FASTME_TYPE=0
FORCE_REBUILD=false

# ── Colors (auto-detect TTY; suppress with NO_COLOR=1) ──────────────────────
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
    C_RESET=$'\033[0m'
    C_BOLD=$'\033[1m'
    C_BLUE=$'\033[34m'
    C_GREEN=$'\033[32m'
    C_YELLOW=$'\033[33m'
    C_RED=$'\033[31m'
    C_CYAN=$'\033[36m'
else
    C_RESET="" ; C_BOLD="" ; C_BLUE="" ; C_GREEN="" ; C_YELLOW="" ; C_RED="" ; C_CYAN=""
fi

log_info()  { echo "${C_BLUE}[STEQ]${C_RESET} $*"; }
log_ok()    { echo "${C_GREEN}[STEQ ✓]${C_RESET} $*"; }
log_warn()  { echo "${C_YELLOW}[STEQ ⚠]${C_RESET} $*" >&2; }
log_err()   { echo "${C_RED}[STEQ ✗] Error:${C_RESET} $*" >&2; }

# ── Usage / Help ─────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
${C_BOLD}STEQ — Species Tree Estimation using Quartet${C_RESET}

${C_BOLD}USAGE${C_RESET}
  $(basename "$0") -i <gene_trees> -o <species_tree> [OPTIONS]

${C_BOLD}DESCRIPTION${C_RESET}
  Infers a species tree from a set of gene trees in Newick format using
  quartet-based distance estimation and a distance-based tree-building method.

${C_BOLD}REQUIRED${C_RESET}
  -i, --input   FILE    Input gene trees file (Newick, one tree per line)
  -o, --output  FILE    Output species tree file

${C_BOLD}OPTIONS${C_RESET}
  -c, --central-tendency  TYPE   Central tendency for distance aggregation
                                   0 = Mean            (default)
                                   1 = Median
                                   2 = Mode
                                   3 = Trimmed Mean    (Normal Distribution, ±2σ)

  -f, --fastme-type       TYPE   Distance-based tree method
                                   0 = FastME_BAL      (balanced minimum evolution, default)
                                   1 = FastME_NJ       (neighbor joining via FastME)
                                   2 = FastME_UNJ      (unweighted NJ via FastME)
                                   Note: if the distance matrix is incomplete,
                                   BIONJ is used automatically regardless of this flag.

  -m, --matrix   FILE    Path for intermediate distance matrix file
                         (default: a temp file is created and cleaned up)

  --rebuild               Force recompilation of STEQ even if binary exists

  -h, --help              Show this help message and exit

${C_BOLD}EXAMPLES${C_RESET}
  # Basic run with default settings (Mean + FastME_BAL)
  ./run_steq.sh -i all_gt_37.tre -o output_37.tre

  # Use Median distances and Neighbor Joining
  ./run_steq.sh -i all_gt_37.tre -o output_37.tre -c 1 -f 1

  # Force rebuild before running
  ./run_steq.sh -i all_gt_37.tre -o output_37.tre --rebuild
EOF
    exit "${1:-0}"
}

# ── Parse Arguments ──────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage 1

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT="$2"; shift 2 ;;
        -o|--output)
            OUTPUT="$2"; shift 2 ;;
        -c|--central-tendency)
            CENTRAL_TENDENCY="$2"; shift 2 ;;
        -f|--fastme-type)
            FASTME_TYPE="$2"; shift 2 ;;
        -m|--matrix)
            MATRIX="$2"; shift 2 ;;
        --rebuild)
            FORCE_REBUILD=true; shift ;;
        -h|--help)
            usage 0 ;;
        *)
            log_err "Unknown option '$1'"
            echo "Run '$(basename "$0") --help' for usage." >&2
            exit 1 ;;
    esac
done

# ── Validate Required Arguments ──────────────────────────────────────────────
if [[ -z "$INPUT" ]]; then
    log_err "Input gene trees file is required (-i)."
    exit 1
fi

if [[ -z "$OUTPUT" ]]; then
    log_err "Output species tree file is required (-o)."
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    log_err "Input file not found: $INPUT"
    exit 1
fi

if [[ "$CENTRAL_TENDENCY" != [0-3] ]]; then
    log_err "Central tendency must be 0, 1, 2, or 3 (got '$CENTRAL_TENDENCY')."
    exit 1
fi

if [[ "$FASTME_TYPE" != [0-2] ]]; then
    log_err "FastME type must be 0, 1, or 2 (got '$FASTME_TYPE')."
    exit 1
fi

# ── Resolve Paths ────────────────────────────────────────────────────────────
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"

OUTPUT_DIR="$(dirname "$OUTPUT")"
mkdir -p "$OUTPUT_DIR" 2>/dev/null || true
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
OUTPUT="$OUTPUT_DIR/$(basename "$OUTPUT")"

CLEANUP_MATRIX=false
if [[ -z "$MATRIX" ]]; then
    MATRIX="$(mktemp /tmp/steq_matrix_XXXXXX.txt)"
    CLEANUP_MATRIX=true
else
    MATRIX_DIR="$(dirname "$MATRIX")"
    mkdir -p "$MATRIX_DIR" 2>/dev/null || true
    MATRIX_DIR="$(cd "$MATRIX_DIR" && pwd)"
    MATRIX="$MATRIX_DIR/$(basename "$MATRIX")"
fi

# ── Check Build Dependencies ────────────────────────────────────────────────
check_dep() {
    if ! command -v "$1" &>/dev/null; then
        log_err "'$1' is required but not found."
        echo "  Install it with:  sudo apt-get install -y $2" >&2
        exit 1
    fi
}

check_dep flex flex
check_dep g++ g++

# ── Build STEQ if Needed ────────────────────────────────────────────────────
if [[ "$FORCE_REBUILD" == true ]] || [[ ! -x "$STEQ_EXEC" ]]; then
    log_info "Compiling STEQ ..."
    pushd "$SUMMARIZER_DIR" > /dev/null
    flex -o _steq_temp.c STEQ.l
    g++ _steq_temp.c -lfl -o STEQ.out
    rm -f _steq_temp.c
    popd > /dev/null
    log_ok "Build complete → $STEQ_EXEC"
else
    log_info "Using existing binary → $STEQ_EXEC"
fi

# Ensure tree-inference binaries are executable
chmod +x "$FASTME_EXEC" "$BIONJ_EXEC" 2>/dev/null || true

# ── Human-Readable Labels ───────────────────────────────────────────────────
CT_LABELS=("Mean" "Median" "Mode" "Trimmed Mean (Normal Distribution)")
FM_LABELS=("FastME_BAL (Balanced Minimum Evolution)" "FastME_NJ (Neighbor Joining)" "FastME_UNJ (Unweighted NJ)")
NUM_TREES="$(wc -l < "$INPUT")"

# ── Print Run Configuration ─────────────────────────────────────────────────
echo ""
echo "${C_BOLD}${C_CYAN}╔══════════════════════════════════════════════════════╗${C_RESET}"
echo "${C_BOLD}${C_CYAN}║           STEQ — Run Configuration                  ║${C_RESET}"
echo "${C_BOLD}${C_CYAN}╚══════════════════════════════════════════════════════╝${C_RESET}"
echo ""
echo "  ${C_BOLD}Input gene trees${C_RESET}   : $INPUT"
echo "  ${C_BOLD}Number of trees${C_RESET}    : $NUM_TREES"
echo "  ${C_BOLD}Central tendency${C_RESET}   : ${CT_LABELS[$CENTRAL_TENDENCY]} ($CENTRAL_TENDENCY)"
echo "  ${C_BOLD}Tree method${C_RESET}        : ${FM_LABELS[$FASTME_TYPE]} ($FASTME_TYPE)"
echo "  ${C_BOLD}Output tree${C_RESET}        : $OUTPUT"
if [[ "$CLEANUP_MATRIX" == false ]]; then
echo "  ${C_BOLD}Distance matrix${C_RESET}    : $MATRIX"
fi
echo ""

# ── Run STEQ ─────────────────────────────────────────────────────────────────
log_info "Running STEQ ..."

# STEQ.out uses relative path ../Binaries/ so we must run from Summarizer/
pushd "$SUMMARIZER_DIR" > /dev/null

START_SEC=$(date +%s)
START_NS="$(date +%s%N 2>/dev/null || echo 0)"

# Suppress harmless FastME stat-file cleanup message while keeping real errors
./STEQ.out "$INPUT" "$CENTRAL_TENDENCY" "$MATRIX" "$OUTPUT" "$FASTME_TYPE" \
    2> >(grep -v 'fastme_stat' >&2)

END_NS="$(date +%s%N 2>/dev/null || echo 0)"
END_SEC=$(date +%s)

popd > /dev/null

# ── Clean Up Temp Matrix ─────────────────────────────────────────────────────
if [[ "$CLEANUP_MATRIX" == true ]]; then
    rm -f "$MATRIX" 2>/dev/null || true
fi

# ── Report Results ───────────────────────────────────────────────────────────
echo ""
if [[ -f "$OUTPUT" ]]; then
    # Compute elapsed time
    if [[ "$START_NS" != "0" && "$END_NS" != "0" ]]; then
        ELAPSED="$(awk -v s="$START_NS" -v e="$END_NS" 'BEGIN { printf "%.3f", (e-s)/1000000000 }')"
    else
        ELAPSED="$((END_SEC - START_SEC))"
    fi

    log_ok "Species tree written to: $OUTPUT"
    echo ""
    echo "  ${C_BOLD}Elapsed time${C_RESET}  : ${ELAPSED}s"
    echo ""
    echo "${C_BOLD}  Output tree (first 200 chars):${C_RESET}"
    echo "  $(head -c 200 "$OUTPUT")"
    echo ""
else
    log_err "Output tree file was not created — check input file format."
    exit 1
fi

