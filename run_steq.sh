#!/usr/bin/env bash
# ==============================================================================
# run_steq.sh — All-inclusive wrapper to build and run STEQ
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUMMARIZER_DIR="$SCRIPT_DIR/Summarizer"
BINARIES_DIR="$SCRIPT_DIR/Binaries"
STEQ_EXEC="$SUMMARIZER_DIR/STEQ.out"
FASTME_EXEC="$BINARIES_DIR/fastme-2.1.5-linux64"
BIONJ_EXEC="$BINARIES_DIR/BIONJ"

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT=""
OUTPUT=""
MATRIX=""
CENTRAL_TENDENCY=0
FASTME_TYPE=0
FORCE_REBUILD=false
NQD_METHOD="baseline"
VALIDATE_NQD=false
BENCHMARK_NQD=false
VALIDATE_PAIRS=20000

# ── Colors (TTY only, can disable with NO_COLOR=1) ──────────────────────────
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
    C_RESET=$'\033[0m'
    C_BOLD=$'\033[1m'
    C_BLUE=$'\033[34m'
    C_GREEN=$'\033[32m'
    C_YELLOW=$'\033[33m'
    C_RED=$'\033[31m'
else
    C_RESET=""
    C_BOLD=""
    C_BLUE=""
    C_GREEN=""
    C_YELLOW=""
    C_RED=""
fi

log_info() { echo "${C_BLUE}[STEQ]${C_RESET} $*"; }
log_ok() { echo "${C_GREEN}[STEQ]${C_RESET} $*"; }
log_warn() { echo "${C_YELLOW}[STEQ]${C_RESET} $*" >&2; }
log_err() { echo "${C_RED}Error:${C_RESET} $*" >&2; }

# ── Usage / Help ──────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: $(basename "$0") -i <gene_trees> -o <species_tree> [OPTIONS]

Run the STEQ (Species Tree Estimation using Quartet) method on a set of
gene trees in Newick format to infer a species tree.

Required:
  -i, --input   FILE    Input gene trees file (Newick, one tree per line)
  -o, --output  FILE    Output species tree file

Options:
  -c, --central-tendency  TYPE   Central tendency measure (default: 0)
                                   0 = Mean
                                   1 = Median
                                   2 = Mode
                                   3 = Normal Distribution (trimmed mean)

  -f, --fastme-type       TYPE   Distance-based tree method (default: 0)
                                   0 = FastME_BAL  (balanced minimum evolution)
                                   1 = FastME_NJ   (neighbor joining via FastME)
                                   2 = FastME_UNJ  (unweighted NJ via FastME)
                                   If the distance matrix is incomplete, BIONJ
                                   is used automatically regardless of this flag.

  -m, --matrix   FILE    Path for intermediate distance matrix file
                         (by default a temp file is used; note: the STEQ
                          binary deletes this file after tree inference)

  --nqd-method           NQD query strategy (default: baseline)
                           baseline  = original path traversal
                           optimized = preprocessed LCA/prefix-sum method

  --validate-nqd         Compare baseline vs optimized NQD values/matrix and
                         stop on mismatch

  --benchmark-nqd        Run baseline and optimized timing benchmark before
                         final tree inference

  --validate-pairs N     Number of random per-tree pairs used in validation
                         when full-pair check is expensive (default: 20000)

  --rebuild               Force recompilation of STEQ even if binary exists

  -h, --help              Show this help message and exit

Examples:
  # Basic run (Mean + FastME_BAL)
  ./run_steq.sh -i gene_trees.tre -o species_tree.tre

  # Median distances, NJ tree building
  ./run_steq.sh -i gene_trees.tre -o species_tree.tre -c 1 -f 1


  # Force recompile
  ./run_steq.sh -i gene_trees.tre -o species_tree.tre --rebuild
EOF
    exit "${1:-0}"
}

# ── Parse Arguments ───────────────────────────────────────────────────────────
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
        --nqd-method)
            NQD_METHOD="$2"; shift 2 ;;
        --validate-nqd)
            VALIDATE_NQD=true; shift ;;
        --benchmark-nqd)
            BENCHMARK_NQD=true; shift ;;
        --validate-pairs)
            VALIDATE_PAIRS="$2"; shift 2 ;;

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

# ── Validate Required Args ───────────────────────────────────────────────────
if [[ -z "$INPUT" ]]; then
    log_err "Input gene trees file is required (-i)."
    exit 1
fi

if [[ -z "$OUTPUT" ]]; then
    log_err "Output species tree file is required (-o)."
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    log_err "Input file '$INPUT' not found."
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

if [[ "$NQD_METHOD" != "baseline" && "$NQD_METHOD" != "optimized" ]]; then
    log_err "--nqd-method must be baseline or optimized (got '$NQD_METHOD')."
    exit 1
fi

if ! [[ "$VALIDATE_PAIRS" =~ ^[0-9]+$ ]]; then
    log_err "--validate-pairs must be a non-negative integer (got '$VALIDATE_PAIRS')."
    exit 1
fi

# ── Resolve Paths (STEQ internals use relative paths from Summarizer/) ────────
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
OUTPUT_DIR="$(cd "$(dirname "$OUTPUT")" 2>/dev/null && pwd)" || mkdir -p "$(dirname "$OUTPUT")"
OUTPUT_DIR="$(cd "$(dirname "$OUTPUT")" && pwd)"
OUTPUT="$OUTPUT_DIR/$(basename "$OUTPUT")"

if [[ -z "$MATRIX" ]]; then
    MATRIX="$(mktemp /tmp/steq_matrix_XXXXXX.txt)"
else
    MATRIX_DIR="$(cd "$(dirname "$MATRIX")" 2>/dev/null && pwd)" || mkdir -p "$(dirname "$MATRIX")"
    MATRIX_DIR="$(cd "$(dirname "$MATRIX")" && pwd)"
    MATRIX="$MATRIX_DIR/$(basename "$MATRIX")"
fi

# ── Check Dependencies ───────────────────────────────────────────────────────
check_dep() {
    if ! command -v "$1" &>/dev/null; then
        echo "Error: '$1' is required but not found. Install it and try again." >&2
        echo "  e.g.  sudo apt-get install -y $2" >&2
        exit 1
    fi
}

check_dep flex flex
check_dep g++ g++

# ── Build STEQ if needed ─────────────────────────────────────────────────────
if [[ "$FORCE_REBUILD" == true ]] || [[ ! -x "$STEQ_EXEC" ]]; then
    log_info "Compiling STEQ.l ..."
    pushd "$SUMMARIZER_DIR" > /dev/null
    flex -o _steq_temp.c STEQ.l
    g++ _steq_temp.c -lfl -o STEQ.out
    rm -f _steq_temp.c
    popd > /dev/null
    log_ok "Build complete: $STEQ_EXEC"
else
    log_info "Using existing binary: $STEQ_EXEC"
fi

# Ensure tree-inference binaries are executable
chmod +x "$FASTME_EXEC" "$BIONJ_EXEC" 2>/dev/null || true

# ── Print Run Summary ────────────────────────────────────────────────────────
CT_LABELS=("Mean" "Median" "Mode" "Normal Distribution")
FM_LABELS=("FastME_BAL" "FastME_NJ" "FastME_UNJ")
NUM_TREES="$(wc -l < "$INPUT")"

echo ""
echo "${C_BOLD}${C_BLUE}[STEQ] ─── Run Configuration ───${C_RESET}"
echo "  Input gene trees : $INPUT  ($NUM_TREES trees)"
echo "  Central tendency : ${CT_LABELS[$CENTRAL_TENDENCY]} ($CENTRAL_TENDENCY)"
echo "  Tree method      : ${FM_LABELS[$FASTME_TYPE]} ($FASTME_TYPE)"
echo "  NQD method       : $NQD_METHOD"
echo "  Validate NQD     : $VALIDATE_NQD"
echo "  Benchmark NQD    : $BENCHMARK_NQD"
echo "  Output tree      : $OUTPUT"
echo ""

# ── Run STEQ ─────────────────────────────────────────────────────────────────
log_info "Running ..."

# STEQ.out uses relative path ../Binaries/ so we must run from Summarizer/
pushd "$SUMMARIZER_DIR" > /dev/null
STEQ_ARGS=(
  "$INPUT"
  "$CENTRAL_TENDENCY"
  "$MATRIX"
  "$OUTPUT"
  "$FASTME_TYPE"
  "--nqd_method" "$NQD_METHOD"
  "--validate_pairs" "$VALIDATE_PAIRS"
)

if [[ "$VALIDATE_NQD" == true ]]; then
    STEQ_ARGS+=("--validate_nqd")
fi

if [[ "$BENCHMARK_NQD" == true ]]; then
    STEQ_ARGS+=("--benchmark_nqd")
fi

./STEQ.out "${STEQ_ARGS[@]}"
popd > /dev/null

# ── Done ──────────────────────────────────────────────────────────────────────
if [[ -f "$OUTPUT" ]]; then
    log_ok "Done! Species tree written to: $OUTPUT"
    echo ""
    echo "${C_BOLD}  Output tree (first 200 chars):${C_RESET}"
    echo "  $(head -c 200 "$OUTPUT")..."
    echo ""
else
    log_warn "Output tree file was not created. Check input file format."
    exit 1
fi
