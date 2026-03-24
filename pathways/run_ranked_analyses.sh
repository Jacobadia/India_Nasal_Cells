#!/usr/bin/env bash

set -euo pipefail

# Path to the ranked file (edit as needed)
RANKED_FILE="../pathway_artifacts/lpm_protein_control_nothing.rnk"

# Path to the GSEA ranked analysis script
GSEA_RANKED_SCRIPT="./gsea_analysis_ranked.sh"

# Directory containing GMT files
GMT_DIR="../pathway_artifacts/msigdb"

# List of GMT files (edit if you want to subset)
GMT_FILES=(
	"h.all.v2026.1.Hs.symbols.gmt"
	"c7.immunesigdb.v2026.1.Hs.symbols.gmt"
	"c7.vax.v2026.1.Hs.symbols.gmt"
	"c2.cp.reactome.v2026.1.Hs.symbols.gmt"
	"c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt"
	"c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"
)

# Output directories for each analysis (edit as needed)
OUTDIRS=(
	"../pathway_artifacts/Hallmark_ranked_GSEA"
	"../pathway_artifacts/ImmuneSigDB_ranked_GSEA"
	"../pathway_artifacts/Vax_ranked_GSEA"
	"../pathway_artifacts/Reactome_ranked_GSEA"
	"../pathway_artifacts/KEGG_Medicus_ranked_GSEA"
	"../pathway_artifacts/KEGG_Legacy_ranked_GSEA"
)

for i in "${!GMT_FILES[@]}"; do
	GMT_FILE="$GMT_DIR/${GMT_FILES[$i]}"
	OUTDIR="${OUTDIRS[$i]}"
	echo "\n==== Running GSEA Preranked for: ${GMT_FILES[$i]} ===="
	bash "$GSEA_RANKED_SCRIPT" --ranked "$RANKED_FILE" --geneset "$GMT_FILE" --outdir "$OUTDIR"
done

echo "\nAll GSEA Preranked analyses complete."
