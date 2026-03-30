#!/usr/bin/env bash

set -euo pipefail

# -----------------------------
# GSEA PRERANKED RUN CONFIGURATION
# Edit these values before each run.
# -----------------------------

GSEA_EXEC_PATH="/mnt/c/Users/bryan/my_scripts/GSEA_4.4.0/gsea-cli.sh"
RANKED_FILE=""
GENESET_FILE=""
OUT_DIR=""

usage() {
	echo "Usage: $0 --ranked <RANKED_FILE.rnk> --geneset <GENESET_FILE.gmt> --outdir <OUT_DIR>"
	echo ""
	echo "Required flags:"
	echo "  -r, --ranked    Path to ranked .rnk file (tab-delimited: gene\tscore)"
	echo "  -g, --geneset   Path to gene set GMT file"
	echo "  -o, --outdir    Output directory for GSEA results"
	echo "  -h, --help      Show this help message"
	echo ""
	echo "Example:"
	echo "  $0 --ranked ../pathway_artifacts/my_ranked.rnk --geneset ../pathway_artifacts/msigdb/h.all.v2026.1.Hs.symbols.gmt --outdir ../pathway_artifacts/Hallmark_GSEA_PRERANKED"
}

while [[ "$#" -gt 0 ]]; do
	case "$1" in
		-r|--ranked)
			RANKED_FILE="${2:-}"
			shift 2
			;;
		-g|--geneset)
			GENESET_FILE="${2:-}"
			shift 2
			;;
		-o|--outdir)
			OUT_DIR="${2:-}"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		*)
			echo "Unknown argument: $1"
			usage
			exit 1
			;;
	esac
done

if [[ -z "$RANKED_FILE" || -z "$GENESET_FILE" || -z "$OUT_DIR" ]]; then
	echo "--ranked, --geneset, and --outdir are all required."
	usage
	exit 1
fi

RPT_LABEL="gsea_preranked_run"
N_PERM="1000"
SEED="149"

if [[ ! -x "$GSEA_EXEC_PATH" ]]; then
	echo "GSEA executable not found or not executable: $GSEA_EXEC_PATH"
	exit 1
fi

if [[ ! -f "$RANKED_FILE" ]]; then
	echo "Ranked file not found: $RANKED_FILE"
	exit 1
fi

if [[ ! -f "$GENESET_FILE" ]]; then
	echo "Gene set file not found: $GENESET_FILE"
	exit 1
fi

mkdir -p "$OUT_DIR"

echo "Running GSEA Preranked with:"
echo "  GSEA CLI : $GSEA_EXEC_PATH"
echo "  RNK      : $RANKED_FILE"
echo "  GMT      : $GENESET_FILE"
echo "  OUT      : $OUT_DIR"
echo "  LABEL    : $RPT_LABEL"
echo "  NPERM    : $N_PERM"
echo "  SEED     : $SEED"

"$GSEA_EXEC_PATH" GSEAPreranked \
	-rnk "$RANKED_FILE" \
	-gmx "$GENESET_FILE" \
	-out "$OUT_DIR" \
	-rpt_label "$RPT_LABEL" \
	-nperm "$N_PERM" \
	-scoring_scheme weighted \
	-set_max 500 \
	-set_min 15 \
	-plot_top_x 20 \
	-rnd_seed "$SEED" \
	-create_svgs false \
	-make_sets true \
	-zip_report false

echo "GSEA Preranked run finished. Results in: $OUT_DIR"
