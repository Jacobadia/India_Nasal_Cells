"""Convert corrected gene counts TSV into a GSEA-ready GCT file.

This script reads the count matrix in ../artifacts/gene_counts_corrected.tsv and
writes a GCT v1.2 file to ../pathway_artifacts without modifying the input file.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from name_mapper import GeneNameMapper


DEFAULT_INPUT = Path("../artifacts/gene_counts_corrected.tsv")
DEFAULT_OUTPUT = Path("../pathway_artifacts/gene_counts_corrected.gct")
MAPPER = GeneNameMapper()

def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Convert a tab-delimited gene count matrix into GSEA GCT format."
	)
	parser.add_argument(
		"--input",
		type=Path,
		default=DEFAULT_INPUT,
		help="Input TSV file (default: ../artifacts/gene_counts_corrected.tsv)",
	)
	parser.add_argument(
		"--output",
		type=Path,
		default=DEFAULT_OUTPUT,
		help="Output GCT file (default: ../pathway_artifacts/gene_counts_corrected.gct)",
	)
	return parser.parse_args()


def convert_tsv_to_gct(input_path: Path, output_path: Path) -> None:
	if not input_path.exists():
		raise FileNotFoundError(f"Input file not found: {input_path}")

	output_path.parent.mkdir(parents=True, exist_ok=True)

	with input_path.open("r", newline="", encoding="utf-8") as infile:
		reader = csv.reader(infile, delimiter="\t")
		header = next(reader, None)
		if not header:
			raise ValueError(f"Input file is empty: {input_path}")

		if len(header) < 7:
			raise ValueError(
				"Input does not have expected columns. "
				"Need at least metadata + one sample column."
			)

		# The source matrix begins with: Geneid, Chr, Start, End, Strand, Length,
		# then one column per sample. GCT needs Name, Description, then sample values.
		sample_names = header[6:]
		rows: list[list[str]] = []

		for row in reader:
			if not row:
				continue

			# Keep only Name (gene id) and expression columns.
			name = row[0] if len(row) > 0 else ""
			expression_values = row[6:] if len(row) > 6 else []
			name = MAPPER.map_gene_id(name)

			if len(expression_values) < len(sample_names):
				expression_values.extend([""] * (len(sample_names) - len(expression_values)))
			elif len(expression_values) > len(sample_names):
				expression_values = expression_values[: len(sample_names)]

			rows.append([name, "na", *expression_values])

	with output_path.open("w", newline="", encoding="utf-8") as outfile:
		writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
		writer.writerow(["#1.2"])
		writer.writerow([str(len(rows)), str(len(sample_names))])
		writer.writerow(["Name", "Description", *sample_names])
		writer.writerows(rows)


def main() -> None:
	args = parse_args()
	convert_tsv_to_gct(args.input, args.output)
	print(f"Wrote GCT file to: {args.output}")


if __name__ == "__main__":
	main()
