#!/usr/bin/env bash

set -euo pipefail

# -----------------------------
# GSEA RUN CONFIGURATION
# Edit these values before each run.
# -----------------------------

GSEA_EXEC_PATH="/mnt/c/Users/bryan/my_scripts/GSEA_4.4.0/gsea-cli.sh"
GCT_FILE="../pathway_artifacts/counts_sex_corrected.gct"
CLS_FILE="../pathway_artifacts/phenotypes.cls"
GENESET_FILE="../pathway_artifacts/msigdb/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt"

# If you need phenotype subset syntax, use e.g. "./phenotypes.cls#Control_vs_Case".
# Leave empty to use CLS_FILE directly.
CLS_ARG=""

# Edit this outdirectory to be what you want.
OUT_DIR="../pathway_artifacts/KEGG_Medicus_GSEA"
RPT_LABEL="gsea_run"
N_PERM="1000"
SEED="149"

if [[ -z "$CLS_ARG" ]]; then
	CLS_ARG="$CLS_FILE"
fi

if [[ -z "$GCT_FILE" || -z "$CLS_FILE" || -z "$GENESET_FILE" ]]; then
	echo "Set GCT_FILE, CLS_FILE, and GENESET_FILE at the top of this script before running."
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
echo "  CLS ARG  : $CLS_ARG"
echo "  GMT      : $GENESET_FILE"
echo "  OUT      : $OUT_DIR"
echo "  LABEL    : $RPT_LABEL"
echo "  NPERM    : $N_PERM"
echo "  SEED     : $SEED"

"$GSEA_EXEC_PATH" GSEA \
	-res "$GCT_FILE" \
	-cls "$CLS_ARG" \
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