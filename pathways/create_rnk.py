from pathlib import Path
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).parent))
from name_mapper import GeneNameMapper

INPUT_PATH = "../artifacts/lpm_protein_control_sex/deg_results_full.txt"
OUTPUT_PATH = "../pathway_artifacts/lpm_protein_control_sex.rnk"


def map_and_filter_genes(df, mapper):
	mapped = []
	not_mapped = 0
	not_protein_coding = 0
	ens_id_name = 0
	for gene_id, row in df.iterrows():
		try:
			gene_name = mapper.map_gene_id(gene_id)
		except Exception:
			not_mapped += 1
			continue
		gene_type = mapper.map_gene_type(gene_id)
		if gene_name.startswith("ENSG"):
			ens_id_name += 1
			continue
		if gene_type != "protein_coding":
			not_protein_coding += 1
			continue
		mapped.append((gene_name, row['t']))
	return mapped, not_mapped, ens_id_name, not_protein_coding

def log_exclusions(not_mapped, ens_id_name, not_protein_coding):
	print(f"Excluded {not_mapped} genes because they could not be mapped.")
	print(f"Excluded {ens_id_name} genes because mapping was to another ENSG ID, not a gene symbol.")
	print(f"Excluded {not_protein_coding} genes because they are not protein_coding.")

def log_duplicates(rnk_df):
	duplicated = rnk_df[rnk_df.duplicated('gene', keep=False)]
	if not duplicated.empty:
		print("Gene names with duplications (and excluded rows):")
		for gene, group in duplicated.groupby('gene'):
			idx_keep = group['t'].abs().idxmax()
			for idx, row in group.iterrows():
				if idx != idx_keep:
					print(f"  {gene}: excluded t={row['t']}")

def deduplicate_genes(rnk_df):
	before_dedup = len(rnk_df)
	log_duplicates(rnk_df)
	rnk_df = rnk_df.loc[rnk_df.groupby('gene')['t'].apply(lambda x: x.abs().idxmax())]
	after_dedup = len(rnk_df)
	print(f"Excluded {before_dedup - after_dedup} genes due to duplicate gene names (kept largest abs(t)).")
	return rnk_df

def main(input_file, output_file):
	df = pd.read_csv(input_file, sep='\t', index_col=0)
	if 't' not in df.columns:
		raise ValueError("Input file must contain a 't' column.")
	mapper = GeneNameMapper()
	mapped, not_mapped, ens_id_name, not_protein_coding = map_and_filter_genes(df, mapper)
	log_exclusions(not_mapped, ens_id_name, not_protein_coding)
	rnk_df = pd.DataFrame(mapped, columns=['gene', 't'])
	rnk_df = deduplicate_genes(rnk_df)
	rnk_df = rnk_df.sort_values('t', ascending=False)
	rnk_df.to_csv(output_file, sep='\t', header=False, index=False)

if __name__ == "__main__":
	main(INPUT_PATH, OUTPUT_PATH)