"""Convert corrected gene counts TSV into a GSEA-ready GCT file.

This script reads the count matrix in ../pathway_artifacts/counts_sex_corrected.tsv
and writes a GCT v1.2 file to ../pathway_artifacts without modifying the input file.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
from name_mapper import GeneNameMapper


# Edit these values before each run.
INPUT_PATH = Path("../pathway_artifacts/counts_sex_corrected.tsv")
OUTPUT_PATH = Path("../pathway_artifacts/counts_sex_corrected.gct")
TARGET_GENE_TYPE = "protein_coding"
MAPPER = GeneNameMapper()


@dataclass
class ConversionStats:
	total_rows: int = 0
	kept_rows: int = 0
	skipped_missing_mapping: int = 0
	skipped_non_protein_coding: int = 0
	skipped_unmapped_symbol: int = 0
	duplicate_gene_symbols: int = 0
	duplicate_gene_symbol_names: list[str] = field(default_factory=list)
	removed_low_variance_duplicates: int = 0


@dataclass
class GctRowVariance:
	"""A GCT row with its calculated variance across samples."""
	row: list[str]
	variance: float


def normalize_expression_values(values: list[str], expected_len: int) -> list[str]:
	if len(values) < expected_len:
		values.extend([""] * (expected_len - len(values)))
	elif len(values) > expected_len:
		values = values[:expected_len]
	return values


def _calculate_row_variance(expression_values: list[str]) -> float:
	"""Calculate variance of expression values across samples."""
	try:
		numeric_values = [float(v) if v else 0.0 for v in expression_values]
	except ValueError:
		numeric_values = [0.0] * len(expression_values)
	
	return float(np.var(numeric_values)) if numeric_values else 0.0


def _process_row(
	row: list[str],
	sample_names: list[str],
	mapper: GeneNameMapper,
	stats: ConversionStats,
) -> tuple[GctRowVariance, str] | None:
	"""Process a single row from the input file.
	
	Returns (GctRowVariance, gene_symbol) or None if row should be skipped.
	"""
	stats.total_rows += 1
	gene_id = row[0].strip() if row and row[0] else ""
	if not gene_id:
		stats.skipped_missing_mapping += 1
		return None

	try:
		gene_type = mapper.map_gene_type(gene_id)
		gene_symbol = mapper.map_gene_id(gene_id).strip()
	except ValueError:
		stats.skipped_missing_mapping += 1
		return None

	if gene_type != TARGET_GENE_TYPE:
		stats.skipped_non_protein_coding += 1
		return None

	if not gene_symbol or gene_symbol.startswith("ENSG"):
		stats.skipped_unmapped_symbol += 1
		return None

	expression_values = normalize_expression_values(row[1:], len(sample_names))
	variance = _calculate_row_variance(expression_values)
	
	gct_row = [gene_symbol, "na", *expression_values]
	return GctRowVariance(gct_row, variance), gene_symbol


def _deduplicate_rows_by_variance(
	rows_by_symbol: defaultdict[str, list[GctRowVariance]],
	kept_gene_symbols: list[str],
	stats: ConversionStats,
) -> list[list[str]]:
	"""Deduplicate rows by keeping the highest variance entry for each gene symbol."""
	final_rows: list[list[str]] = []
	gene_symbol_counts = Counter(kept_gene_symbols)
	
	for gene_symbol, row_variance_list in rows_by_symbol.items():
		# Keep the row with highest variance
		best_entry = max(row_variance_list, key=lambda entry: entry.variance)
		final_rows.append(best_entry.row)
		
		# Track duplicates and removed rows
		if gene_symbol_counts[gene_symbol] > 1:
			stats.duplicate_gene_symbol_names.append(gene_symbol)
			stats.removed_low_variance_duplicates += gene_symbol_counts[gene_symbol] - 1
	
	stats.duplicate_gene_symbols = len(stats.duplicate_gene_symbol_names)
	stats.duplicate_gene_symbol_names.sort()
	stats.kept_rows = len(final_rows)
	
	return final_rows


def build_gct_rows(
	reader: csv.reader,
	sample_names: list[str],
	mapper: GeneNameMapper,
) -> tuple[list[list[str]], ConversionStats]:
	rows_by_symbol: defaultdict[str, list[GctRowVariance]] = defaultdict(list)
	stats = ConversionStats()
	kept_gene_symbols: list[str] = []

	for row in reader:
		if not row:
			continue

		result = _process_row(row, sample_names, mapper, stats)
		if result is None:
			continue
		
		row_variance, gene_symbol = result
		rows_by_symbol[gene_symbol].append(row_variance)
		kept_gene_symbols.append(gene_symbol)
		stats.kept_rows += 1

	final_rows = _deduplicate_rows_by_variance(rows_by_symbol, kept_gene_symbols, stats)
	return final_rows, stats


def convert_tsv_to_gct(input_path: Path, output_path: Path) -> ConversionStats:
	if not input_path.exists():
		raise FileNotFoundError(f"Input file not found: {input_path}")

	output_path.parent.mkdir(parents=True, exist_ok=True)

	with input_path.open("r", newline="", encoding="utf-8") as infile:
		reader = csv.reader(infile, delimiter="\t")
		header = next(reader, None)
		if not header:
			raise ValueError(f"Input file is empty: {input_path}")

		if len(header) < 2:
			raise ValueError(
				"Input does not have expected columns. "
				"Need at least Geneid + one sample column."
			)

		sample_names = header[1:]
		rows, stats = build_gct_rows(reader, sample_names, MAPPER)

	with output_path.open("w", newline="", encoding="utf-8") as outfile:
		writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
		writer.writerow(["#1.2"])
		writer.writerow([str(len(rows)), str(len(sample_names))])
		writer.writerow(["Name", "Description", *sample_names])
		writer.writerows(rows)

	return stats


def main() -> None:
	stats = convert_tsv_to_gct(INPUT_PATH, OUTPUT_PATH)
	print(f"Wrote GCT file to: {OUTPUT_PATH}")
	print(
		"Rows processed: "
		f"{stats.total_rows}\n"
		f"Kept (after deduplication): {stats.kept_rows}\n"
		f"Skipped missing mapping: {stats.skipped_missing_mapping}\n"
		f"Skipped non-protein_coding: {stats.skipped_non_protein_coding}\n"
		f"Skipped ENSG/unmapped symbols: {stats.skipped_unmapped_symbol}\n"
		f"Removed low-variance duplicates: {stats.removed_low_variance_duplicates}\n"
		f"Duplicate kept gene symbols: {stats.duplicate_gene_symbols}"
	)
	if stats.duplicate_gene_symbol_names:
		print("Duplicated symbols (kept highest variance):")
		for symbol in stats.duplicate_gene_symbol_names:
			print(f"  {symbol}")


if __name__ == "__main__":
	main()
