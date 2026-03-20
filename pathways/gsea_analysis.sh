#!/usr/bin/env bash

set -euo pipefail

# -----------------------------
# GSEA RUN CONFIGURATION
# Edit these values before each run.
# -----------------------------

GSEA_EXEC_PATH="/mnt/c/Users/bryan/my_scripts/GSEA_4.4.0/gsea-cli.sh"
GCT_FILE="../pathway_artifacts/counts_sex_corrected.gct"
CLS_FILE="../pathway_artifacts/phenotypes.cls"

GENESET_FILE=""
OUT_DIR=""

usage() {
	echo "Usage: $0 --geneset <GENESET_FILE.gmt> --outdir <OUT_DIR>"
	echo ""
	echo "Required flags:"
	echo "  -g, --geneset   Path to gene set GMT file"
	echo "  -o, --outdir    Output directory for GSEA results"
	echo "  -h, --help      Show this help message"
	echo ""
	echo "Example:"
	echo "  $0 --geneset ../pathway_artifacts/msigdb/h.all.v2026.1.Hs.symbols.gmt --outdir ../pathway_artifacts/Hallmark_GSEA"
}

while [[ "$#" -gt 0 ]]; do
	case "$1" in
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

if [[ -z "$GENESET_FILE" || -z "$OUT_DIR" ]]; then
	echo "Both --geneset and --outdir are required."
	usage
	exit 1
fi

RPT_LABEL="gsea_run"
N_PERM="1000"
SEED="149"

if [[ -z "$GCT_FILE" || -z "$CLS_FILE" || -z "$GENESET_FILE" ]]; then
	echo "Set GCT_FILE and CLS_FILE at the top of this script, and pass --geneset when running."
	exit 1
fi

if [[ ! -x "$GSEA_EXEC_PATH" ]]; then
	echo "GSEA executable not found or not executable: $GSEA_EXEC_PATH"
	exit 1
fi

if [[ ! -f "$GCT_FILE" ]]; then
	echo "GCT file not found: $GCT_FILE"
	exit 1
fi

if [[ ! -f "$CLS_FILE" ]]; then
	echo "CLS file not found: $CLS_FILE"
	exit 1
fi

if [[ ! -f "$GENESET_FILE" ]]; then
	echo "Gene set file not found: $GENESET_FILE"
	exit 1
fi

mkdir -p "$OUT_DIR"

echo "Running GSEA with:"
echo "  GSEA CLI : $GSEA_EXEC_PATH"
echo "  GCT      : $GCT_FILE"
echo "  CLS      : $CLS_FILE"
echo "  GMT      : $GENESET_FILE"
echo "  OUT      : $OUT_DIR"
echo "  LABEL    : $RPT_LABEL"
echo "  NPERM    : $N_PERM"
echo "  SEED     : $SEED"

"$GSEA_EXEC_PATH" GSEA \
	-res "$GCT_FILE" \
	-cls "$CLS_FILE" \
	-gmx "$GENESET_FILE" \
	-out "$OUT_DIR" \
	-rpt_label "$RPT_LABEL" \
	-collapse No_Collapse \
	-mode Max_probe \
	-norm meandiv \
	-nperm "$N_PERM" \
	-permute phenotype \
	-scoring_scheme weighted \
	-set_max 500 \
	-set_min 15 \
	-plot_top_x 20 \
	-rnd_seed "$SEED" \
	-create_svgs false \
	-make_sets true \
	-zip_report false

echo "GSEA run finished. Results in: $OUT_DIR"