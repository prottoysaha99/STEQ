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

        --rebuild)
            FORCE_REBUILD=true; shift ;;
        -h|--help)
            usage 0 ;;
        *)
            echo "Error: Unknown option '$1'" >&2
            echo "Run '$(basename "$0") --help' for usage." >&2
            exit 1 ;;
    esac
done

# ── Validate Required Args ───────────────────────────────────────────────────
if [[ -z "$INPUT" ]]; then
    echo "Error: Input gene trees file is required (-i)." >&2
    exit 1
fi

if [[ -z "$OUTPUT" ]]; then
    echo "Error: Output species tree file is required (-o)." >&2
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file '$INPUT' not found." >&2
    exit 1
fi

if [[ "$CENTRAL_TENDENCY" != [0-3] ]]; then
    echo "Error: Central tendency must be 0, 1, 2, or 3 (got '$CENTRAL_TENDENCY')." >&2
    exit 1
fi

if [[ "$FASTME_TYPE" != [0-2] ]]; then
    echo "Error: FastME type must be 0, 1, or 2 (got '$FASTME_TYPE')." >&2
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
    echo "[STEQ] Compiling STEQ.l ..."
    pushd "$SUMMARIZER_DIR" > /dev/null
    flex -o _steq_temp.c STEQ.l
    g++ _steq_temp.c -lfl -o STEQ.out
    rm -f _steq_temp.c
    popd > /dev/null
    echo "[STEQ] Build complete: $STEQ_EXEC"
else
    echo "[STEQ] Using existing binary: $STEQ_EXEC"
fi

# Ensure tree-inference binaries are executable
chmod +x "$FASTME_EXEC" "$BIONJ_EXEC" 2>/dev/null || true

# ── Print Run Summary ────────────────────────────────────────────────────────
CT_LABELS=("Mean" "Median" "Mode" "Normal Distribution")
FM_LABELS=("FastME_BAL" "FastME_NJ" "FastME_UNJ")
NUM_TREES="$(wc -l < "$INPUT")"

echo ""
echo "[STEQ] ─── Run Configuration ───"
echo "  Input gene trees : $INPUT  ($NUM_TREES trees)"
echo "  Central tendency : ${CT_LABELS[$CENTRAL_TENDENCY]} ($CENTRAL_TENDENCY)"
echo "  Tree method      : ${FM_LABELS[$FASTME_TYPE]} ($FASTME_TYPE)"
echo "  Output tree      : $OUTPUT"
echo ""

# ── Run STEQ ─────────────────────────────────────────────────────────────────
echo "[STEQ] Running ..."

# STEQ.out uses relative path ../Binaries/ so we must run from Summarizer/
pushd "$SUMMARIZER_DIR" > /dev/null
./STEQ.out "$INPUT" "$CENTRAL_TENDENCY" "$MATRIX" "$OUTPUT" "$FASTME_TYPE"
popd > /dev/null

# ── Done ──────────────────────────────────────────────────────────────────────
if [[ -f "$OUTPUT" ]]; then
    echo "[STEQ] Done! Species tree written to: $OUTPUT"
    echo ""
    echo "  Output tree (first 200 chars):"
    echo "  $(head -c 200 "$OUTPUT")..."
    echo ""
else
    echo "[STEQ] Warning: Output tree file was not created. Check input file format." >&2
    exit 1
fi
