#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MSIGDB_DIR="$SCRIPT_DIR/../pathway_artifacts/msigdb"
OUT_BASE_DIR="$SCRIPT_DIR/../pathway_artifacts"
GSEA_SCRIPT="$SCRIPT_DIR/gsea_analysis_nonranked.sh"

if [[ ! -x "$GSEA_SCRIPT" ]]; then
	echo "gsea_analysis_nonranked.sh not found or not executable: $GSEA_SCRIPT"
	exit 1
fi

if [[ ! -d "$MSIGDB_DIR" ]]; then
	echo "MSigDB directory not found: $MSIGDB_DIR"
	exit 1
fi

map_outdir_name() {
	local gmt_name="$1"
	case "$gmt_name" in
		h.all.v*.Hs.symbols.gmt)
			echo "Hallmark_GSEA"
			;;
		c7.immunesigdb.v*.Hs.symbols.gmt)
			echo "ImmuneSigDB_GSEA"
			;;
		c2.cp.kegg_legacy.v*.Hs.symbols.gmt)
			echo "KEGG_Legacy_GSEA"
			;;
		c2.cp.kegg_medicus.v*.Hs.symbols.gmt)
			echo "KEGG_Medicus_GSEA"
			;;
		c2.cp.reactome.v*.Hs.symbols.gmt)
			echo "Reactome_GSEA"
			;;
		c7.vax.v*.Hs.symbols.gmt)
			echo "Vax_GSEA"
			;;
		*)
			local stem="${gmt_name%.gmt}"
			stem="${stem//./_}"
			echo "${stem}_GSEA"
			;;
	esac
}

mapfile -t gmt_files < <(find "$MSIGDB_DIR" -maxdepth 1 -type f -name "*.gmt" | sort)

if [[ "${#gmt_files[@]}" -eq 0 ]]; then
	echo "No .gmt files found in: $MSIGDB_DIR"
	exit 1
fi

echo "Found ${#gmt_files[@]} gene set file(s) in: $MSIGDB_DIR"

for geneset_file in "${gmt_files[@]}"; do
	geneset_name="$(basename "$geneset_file")"
	out_dir_name="$(map_outdir_name "$geneset_name")"
	out_dir="$OUT_BASE_DIR/$out_dir_name"

	echo
	echo "Running GSEA for: $geneset_name"
	echo "Output directory: $out_dir"

	bash "$GSEA_SCRIPT" --geneset "$geneset_file" --outdir "$out_dir"
done

echo
echo "All pathway analyses completed."
