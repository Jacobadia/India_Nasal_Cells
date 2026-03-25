

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


def get_max_name_length(table_data):
    """Return the maximum gene set name length across all tables."""
    return max(df['NAME'].astype(str).map(len).max() for _, df in table_data)

def wrap_gene_set_names(df, wrap_limit=40):
    """Return a copy of df with NAME column unchanged (no wrapping)."""
    return df.copy()


# Helper for table rendering: color for 'Upregulated in' column
def get_upregulated_color(val):
    """Return color for 'Upregulated in' column: red for 'Active', blue for 'Latent', black otherwise."""
    if val == 'Active':
        return '#d62728'  # red
    elif val == 'Latent':
        return '#1f77b4'  # blue
    return 'black'

def render_table_on_axis(ax, df, title, wrap_limit=40):
    ax.axis('off')
    df_disp = wrap_gene_set_names(df, wrap_limit)
    tbl = ax.table(cellText=df_disp.values,
                   colLabels=df_disp.columns,
                   loc='center',
                   cellLoc='center',
                   colLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    # Set each column to its own minimum width and apply cell formatting
    for (row, col), cell in tbl.get_celld().items():
        apply_cell_formatting(cell, row, col, df_disp)
    tbl.scale(1, 1.5)
    ax.set_title(title, fontsize=14, fontweight='bold', pad=10)


def plot_gsea_tables(table_data, output_path=OUTPUT_PATH):
    """Plots all tables in a single matplotlib figure and saves to output_path. Ensures gene set names fit."""
    n_tables = len(table_data)
    n_cols = len(table_data[0][1].columns)
    row_height = 0.5
    max_name_len = get_max_name_length(table_data)
    base_width = 2 + n_cols * 2
    extra_width = min(max((max_name_len - 20) * 0.15, 0), 10)
    fig_width = min(base_width + extra_width, 24)
    fig_height = n_tables * (10 * row_height + 1)
    fig, axes = plt.subplots(n_tables, 1, figsize=(fig_width, fig_height))
    if n_tables == 1:
        axes = [axes]
    for ax, (title, df) in zip(axes, table_data):
        render_table_on_axis(ax, df, title, wrap_limit=40)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved publication-style table as {output_path}")


def main():
    table_data = collect_gsea_tables(GSEA_DIRS)
    plot_gsea_tables(table_data, output_path=OUTPUT_PATH)


#########################
# Table rendering helpers
#########################
def calculate_column_widths(df_disp):
    """Return a list of minimum widths for each column based on the widest cell (header or value)."""
    col_widths = []
    for col_idx, col_name in enumerate(df_disp.columns):
        max_len = max([len(str(col_name))] + [len(str(val)) for val in df_disp.iloc[:, col_idx]])
        # Reduced scale: 0.011 per character, min 0.08, max 0.35
        width = min(max(0.011 * max_len + 0.08, 0.08), 0.35)
        col_widths.append(width)
    return col_widths

def apply_cell_formatting(cell, row, col, df_disp):
    """Set width and color for a table cell."""
    col_widths = calculate_column_widths(df_disp)
    if col < len(col_widths):
        cell.set_width(col_widths[col])
    # Apply color to 'Upregulated in' column
    if row > 0 and col < len(df_disp.columns) and df_disp.columns[col] == 'Upregulated in':
        val = df_disp.iloc[row-1, col]
        color = get_upregulated_color(val)
        cell.set_text_props(color=color)

if __name__ == "__main__":
    main()