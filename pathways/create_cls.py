"""Generate a GSEA-ready CLS file from a GCT matrix.

This script reads sample names from the GCT header and writes a categorical
CLS file with two classes:
- sample names ending in 'LA' -> active
- sample names ending in 'LB' -> latent
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DEFAULT_INPUT = Path("../pathway_artifacts/gene_counts_corrected.gct")
DEFAULT_OUTPUT = Path("../pathway_artifacts/phenotypes.cls")


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Generate a GSEA CLS file from a GCT file."
	)
	parser.add_argument(
		"--input",
		type=Path,
		default=DEFAULT_INPUT,
		help=f"Input GCT file (default: {DEFAULT_INPUT})",
	)
	parser.add_argument(
		"--output",
		type=Path,
		default=DEFAULT_OUTPUT,
		help=f"Output CLS file (default: {DEFAULT_OUTPUT})",
	)
	return parser.parse_args()


def labels_from_sample_names(sample_names: list[str]) -> list[str]:
	labels: list[str] = []
	unrecognized: list[str] = []

	for sample in sample_names:
		clean_name = sample.strip()
		upper_name = clean_name.upper()
		if upper_name.endswith("LA"):
			labels.append("active")
		elif upper_name.endswith("LB"):
			labels.append("latent")
		else:
			unrecognized.append(clean_name)

	if unrecognized:
		joined = ", ".join(unrecognized)
		raise ValueError(
			"Could not assign class labels for sample(s): "
			f"{joined}. Expected names ending in 'LA' or 'LB'."
		)

	return labels


def read_gct_sample_names(gct_path: Path) -> list[str]:
	if not gct_path.exists():
		raise FileNotFoundError(f"Input GCT file not found: {gct_path}")

	with gct_path.open("r", newline="", encoding="utf-8") as infile:
		reader = csv.reader(infile, delimiter="\t")

		version_row = next(reader, None)
		dims_row = next(reader, None)
		header_row = next(reader, None)

		if not version_row or not dims_row or not header_row:
			raise ValueError(f"GCT appears incomplete: {gct_path}")

		if not version_row[0].startswith("#1"):
			raise ValueError(
				f"Unrecognized GCT version in first line: {version_row[0]}"
			)

		if len(header_row) < 3:
			raise ValueError(
				"GCT header is missing sample columns. Expected at least Name, "
				"Description, and one sample column."
			)

		return header_row[2:]


def write_cls(sample_names: list[str], labels: list[str], output_path: Path) -> None:
	output_path.parent.mkdir(parents=True, exist_ok=True)

	with output_path.open("w", encoding="utf-8", newline="\n") as outfile:
		outfile.write(f"{len(sample_names)} 2 1\n")
		outfile.write("# active latent\n")
		outfile.write(" ".join(labels) + "\n")


def generate_cls(input_gct: Path, output_cls: Path) -> None:
	sample_names = read_gct_sample_names(input_gct)
	labels = labels_from_sample_names(sample_names)
	write_cls(sample_names, labels, output_cls)


def main() -> None:
	args = parse_args()
	generate_cls(args.input, args.output)
	print(f"Wrote CLS file to: {args.output}")


if __name__ == "__main__":
	main()