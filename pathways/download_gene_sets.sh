#!/usr/bin/env bash

set -euo pipefail

# -----------------------------
# MSigDB DOWNLOAD CONFIGURATION
# Edit these values before running.
# -----------------------------

MSIGDB_RELEASE="2026.1.Hs"

# Output folder for downloaded GMT files.
OUT_DIR="../pathway_artifacts/msigdb"

# Base URL for downloadable MSigDB release files.
BASE_URL="https://data.broadinstitute.org/gsea-msigdb/msigdb/release"

# Add only the collections you want.
# TODO: Ask Dr. Kahmbati which kegg. There's kegg_medicus and kegg_legacy
COLLECTIONS=(
	"h.all"
	"c7.immunesigdb"
    "c7.vax"
    "c2.cp.reactome"
    "c2.cp.kegg_medicus"
    "c2.cp.kegg_legacy"
)

# Identifier space: symbols or entrez
ID_TYPE="symbols"

mkdir -p "$OUT_DIR"

echo "MSigDB release  : $MSIGDB_RELEASE"
echo "Output directory: $OUT_DIR"
echo "ID type         : $ID_TYPE"

for collection in "${COLLECTIONS[@]}"; do
	filename="${collection}.v${MSIGDB_RELEASE}.${ID_TYPE}.gmt"
	url="${BASE_URL}/${MSIGDB_RELEASE}/${filename}"
	out_file="${OUT_DIR}/${filename}"

	echo "Downloading ${filename}"

	if ! wget -q --show-progress -O "$out_file" "$url"; then
		echo "Failed to download: $url"
		rm -f "$out_file"
		continue
	fi

	if [[ ! -s "$out_file" ]]; then
		echo "Downloaded file is empty: $out_file"
		rm -f "$out_file"
		continue
	fi

	echo "Saved: $out_file"
done

echo "Done."

