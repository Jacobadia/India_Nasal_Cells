# India_Nasal_Cells
This repository contains data and code used in the analysis of the RNA sequencing nasal Indian data

# Code (In Order)

1. **Pre-Trim Quality Control**
Package(s) needed:
- fastqc
- multiqc

Program to run:
- `initial_check.sh`

2. **Trimming and Post-Trim Quality Control**
- `trim_N_check.sh`
3. **STAR Alignment**
- `alignment_alt.sh`
4. **Generate Feature Counts Table**
- `featureCounts.sh`
5. **Limma Analysis**
- Ensure that you are running a POSIX compliant shell (e.g. bash, zsh, etc.) and have R installed on your system. Specifically, Rscript should be on PATH.
- Install the following R packages:
    - CRAN: ggplot2
    - Bioconductor: DESeq2, limma, edgeR

  If you need to install them, use:

  ```r
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
  }

  install.packages("ggplot2")
  BiocManager::install(c("DESeq2", "limma", "edgeR"))
  ```

- In the root project directory, create (if it doesn't exist), a folder named
"artifacts". In this folder, place a copy of the following files: 
    - gene_counts_corrected.tsv (This is the file generated from the STAR alignment and feature counts steps)
    - genetype_lookup.txt (This is the file that we generated from the genome annotatios gtf file. We'll use this to filter out non-protein coding genes from our analysis and map gene ids to gene names)
    - metadata.tsv (This contains the metadata, including the sample names, disease status, and age of the patients. This is used in the limma analysis when we control for sex later).
- cd into the deg folder with:

```bash
cd deg
```

- Run the following code in your shell to perform the Differential Gene Expression analysis:

```bash
Rscript deg/deg_analyses.R
```

- Note: we have commented out all the code for the analyses that we didn't use in our results in the interest of saving you time (It will take another 10-15 minutes if you choose to uncomment the other tests). If you want to run the other analyses too, simply uncomment the code for those analyses found at the very end of the script and run the script again.
- If it succeeded, in the artifacts directory, there should be a lpm_protein_control_nothing/ and an lpm_protein_control_sex/ directory.
- In those, you'll see the pvalue histogram and volcano plot (used in our paper). Additionally, you'll see the full results table order by FDR corrected pvalue. We will use the t statistic from the control sex results table as our ranking metric for the GSEA pathway analysis later on.
- In the lpm_protein_control_nothing/ results significant File, you will find our two primary DEGs, which we used in the subsequent machine learning step.
6. **TBSignatureProfiler**
7. **ROC Pipeline**
8. **GSEA Pathway Analysis**
- Ensure that you are running a POSIX compliant shell (e.g. bash, zsh, etc.). The wget utility should be installed on your system and be on PATH. We will be using it to download the gene sets from the MSigDB database. You can check if you have wget installed by running:

```bash
wget --version
```

- Ensure that Java in installed on your system and that the version is 21 or later. You can check your Java version by running:

```bash
java -version
```

- Ensure that you have the GSEA software (invoked via the CLI) on your system. you can download it from the Msigdb website here: https://www.gsea-msigdb.org/gsea/downloads.jsp. It may prompt you to input your email. Choose the version:
    - GSEA v4.4.0 for the command line (all platforms)
- Unzip the GSEA distibution. It should look something like this.
![IMG](reproduce_utils/gsea_dist.png)
- In the unzipped GSEA distribution, there should be a file called "gsea-cli.sh". This file is VERY IMPORTANT. save the ABSOLUTE PATH to this file. We will use it to invoke the GSEA software from the command line.

- In the pathways/gsea_analysis_ranked.sh, there is a variable at the top called GSEA_EXEC_PATH. Set this variable to the absolute path of the gsea-cli.sh file that we just talked about. It should look something like this:

```bash
GSEA_EXEC_PATH="/path/to/gsea-cli.sh"
```

- Ensure that python is installed an on PATH. Ensure that the following python packages are installed:
    - pandas
    - matplotlib
    - numpy

You can install these packages using pip if you don't have them already:

```bash
pip install pandas matplotlib numpy
```

- Now we can run the analysis.
- cd into the pathways directory with:

```bash
cd pathways
```

- Run the following command to download the gene set files from msigdb.

```bash
./download_gene_sets.sh
```

- There should now be a folder pathway_artifacts/msigdb containing all the gene set files that we will use in our analysis.

- run the following command to run the GSEA analysis:
    - Note, this should take like 10 minutes if your computer is slow like mine

```bash
./run_ranked_analyses.sh
```

- Generate the pathway analysis figure with the following command:

```bash
python create_gsea_figure.py
```

- There should now be a figure called gsea_publication_figure.png in the pathway_artifacts directory. This is the figure that we used in our paper.

9. **Figure Creation**