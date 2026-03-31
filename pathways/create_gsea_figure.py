import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import textwrap

os.chdir(os.path.dirname(os.path.abspath(__file__)))

OUTPUT_PATH = "../pathway_artifacts/gsea_publication_table.png"

GSEA_DIRS = (
    "../pathway_artifacts/Hallmark_ranked_sex_GSEA",
    "../pathway_artifacts/ImmuneSigDB_ranked_sex_GSEA",
    "../pathway_artifacts/KEGG_Legacy_ranked_sex_GSEA",
    "../pathway_artifacts/KEGG_Medicus_ranked_sex_GSEA",
    "../pathway_artifacts/Reactome_ranked_sex_GSEA",
    "../pathway_artifacts/Vax_ranked_sex_GSEA",
)

def find_gsea_reports(gsea_dirs):
    results = []
    for dir in gsea_dirs:
        dir_path = Path(dir)
        assert dir_path.exists(), f"{dir} does not exist"

        # Find the single GSEA result folder inside this directory
        subdirs = [d for d in dir_path.iterdir() if d.is_dir()]
        assert len(subdirs) == 1, f"{dir} has {len(subdirs)} subdirectories (expected 1): {[str(d) for d in subdirs]}"
        result_dir = subdirs[0]

        # Find the .tsv files for both phenotypes
        pos_files = list(result_dir.glob("gsea_report_for_na_pos_*.tsv"))
        neg_files = list(result_dir.glob("gsea_report_for_na_neg_*.tsv"))
        assert len(pos_files) == 1, f"{result_dir} has {len(pos_files)} pos .tsv files (expected 1): {[str(f) for f in pos_files]}"
        assert len(neg_files) == 1, f"{result_dir} has {len(neg_files)} neg .tsv files (expected 1): {[str(f) for f in neg_files]}"

        results.append((dir_path.name, str(pos_files[0]), str(neg_files[0])))
    return results


def get_top10_unique_by_fdr(pos_path, neg_path):
    """
    Given paths to na_pos and na_neg GSEA report TSVs, return a DataFrame of the top 10 unique entries
    (by NAME) with the lowest FDR q-val, considering both files.
    """
    pos_df = pd.read_csv(pos_path, sep='\t')
    neg_df = pd.read_csv(neg_path, sep='\t')
    pos_top = pos_df.nsmallest(10, 'FDR q-val')
    neg_top = neg_df.nsmallest(10, 'FDR q-val')
    combined = pd.concat([pos_top, neg_top], ignore_index=True)
    combined = combined.drop_duplicates(subset=['NAME'])
    combined = combined.nsmallest(10, 'FDR q-val')
    return combined

def collect_gsea_tables(gsea_dirs):
    """Collects and returns a list of (title, DataFrame) for each gene set db, with only selected columns and upregulation info."""
    reports = find_gsea_reports(gsea_dirs)
    table_data = []
    for dbname, pos_path, neg_path in reports:
        top10 = get_top10_unique_by_fdr(pos_path, neg_path)
        top10 = top10.sort_values('FDR q-val', ascending=True)
        # Use original FDR q-val values for full precision
        # Add 'Upregulated in' column
        upreg = top10['NES'].apply(lambda x: 'Active' if x > 0 else 'Latent')
        top10 = top10.assign(**{'Upregulated in': upreg})
        # Only keep NAME, FDR q-val, NES, Upregulated in
        top10 = top10[['NAME', 'FDR q-val', 'NES', 'Upregulated in']]
        table_data.append((dbname, top10))
    return table_data



def main():
    table_data = collect_gsea_tables(GSEA_DIRS)
    print(table_data)

if __name__ == "__main__":
    main()