import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import textwrap

os.chdir(os.path.dirname(os.path.abspath(__file__)))

OUTPUT_PATH = "../pathway_artifacts/gsea_publication_figure.png"

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




import numpy as np
import matplotlib as mpl

def setup_figure(n_panels, ncols, nrows, figsize):
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    plt.subplots_adjust(wspace=0.5, hspace=0.6)  # more vertical space
    return fig, axes, nrows, ncols

def get_global_nes_limits(table_data):
    all_nes = pd.concat([df["NES"] for _, df in table_data])
    nes_min = all_nes.min()
    nes_max = all_nes.max()
    nes_abs = max(abs(nes_min), abs(nes_max))
    return (-nes_abs, nes_abs)

def plot_panel(ax, dbname, df, nes_lim, upreg_palette):
    df_sorted = df.sort_values("NES", ascending=True)
    def snake_to_title(s):
        return s.replace("_", " ").title()
    y_labels = df_sorted["NAME"].apply(snake_to_title)
    y_pos = np.arange(len(y_labels))
    min_padj = 1e-6
    capped_log10 = -np.log10(df_sorted["FDR q-val"].replace(0, min_padj))
    capped_log10 = np.clip(capped_log10, None, 10)
    sizes = capped_log10 * 100
    colors = df_sorted["Upregulated in"].map(upreg_palette)
    scatter = ax.scatter(
        df_sorted["NES"],
        y_pos,
        s=sizes,
        c=colors,
        alpha=0.8,
        edgecolor="k",
        linewidth=0.5
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels([textwrap.fill(n, 60) for n in y_labels])
    ax.set_xlabel("Normalized Enrichment Score (NES)")
    ax.set_title(dbname.replace("_ranked_sex_GSEA", ""))
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    ax.set_xlim(nes_lim)

def get_upregulation_legend_handles(upreg_palette):
    handles = [mpl.lines.Line2D([0], [0], marker='o', color='w', label=lab,
                                 markerfacecolor=col, markersize=10, markeredgecolor='k')
               for lab, col in upreg_palette.items()]
    labels = list(upreg_palette.keys())
    return handles, labels

def get_dot_size_legend_handles():
    padj_legend = [0.05, 0.01, 0.001, 1e-6]
    size_legend = [-np.log10(p) * 100 for p in padj_legend]
    labels = [f"padj={p}" for p in padj_legend[:-1]] + [r"padj < 1e-6"]
    handles = [plt.scatter([], [], s=s, c='gray', alpha=0.7, edgecolor='k') for s in size_legend]
    return handles, labels

def main():
    table_data = collect_gsea_tables(GSEA_DIRS)
    n_panels = len(table_data)
    fig, axes, nrows, ncols = setup_figure(n_panels, ncols=1, nrows=6, figsize=(18, 15))
    upreg_palette = {"Active": "#1f77b4", "Latent": "#d62728"}
    nes_lim = get_global_nes_limits(table_data)

    for idx, (dbname, df) in enumerate(table_data):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        plot_panel(ax, dbname, df, nes_lim, upreg_palette)

    # Remove empty subplots
    for idx in range(n_panels, nrows * ncols):
        row, col = divmod(idx, ncols)
        fig.delaxes(axes[row][col])



    # Separate legends: Upregulation and Dot Size, both at top right, side by side
    upreg_handles, upreg_labels = get_upregulation_legend_handles(upreg_palette)
    dot_handles, dot_labels = get_dot_size_legend_handles()


    # Upregulation legend (top right, shifted further down)
    leg1 = fig.legend(
        handles=upreg_handles,
        labels=[f"Upregulated in: {lab}" for lab in upreg_labels],
        loc="upper right",
        bbox_to_anchor=(0.98, 0.995),  # shifted further down
        frameon=True,
        fontsize=10,
        borderaxespad=0.5,
        bbox_transform=fig.transFigure
    )

    # Dot size legend (vertical, left of upregulation legend, same height)
    leg2 = fig.legend(
        handles=dot_handles,
        labels=[f"Dot size: {lab}" for lab in dot_labels],
        loc="upper right",
        bbox_to_anchor=(0.20, 0.995),  # shifted further down
        frameon=True,
        fontsize=10,
        borderaxespad=0.5,  # increased for a taller legend box
        bbox_transform=fig.transFigure,
        ncol=1,  # vertical legend
        labelspacing=1.0,  # moderate spacing between labels
        handleheight=1.5,  # moderate height for each handle
        handletextpad=1.0  # moderate space between handle and text
    )

    # Make sure both legends are drawn
    fig.add_artist(leg1)
    fig.add_artist(leg2)


    fig.suptitle("GSEA Results: Top Pathways per Database", fontsize=18, y=0.99)
    # Reserve even more space at the top for legends
    plt.tight_layout(rect=[0, 0, 0.97, 0.92])
    plt.savefig(OUTPUT_PATH, dpi=300)
    print(f"Saved GSEA figure to {OUTPUT_PATH}")

if __name__ == "__main__":
    main()