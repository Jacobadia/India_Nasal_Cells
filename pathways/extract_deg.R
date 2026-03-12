source("../methods/map_gene_name.R")

results_file = "../artifacts/pure_control_sex/deg_results_full.txt"
FDR_threshold = 0.15

results = read.delim(results_file)
results_sig = results[results$padj < FDR_threshold, ]

ensemble_ids = rownames(results_sig)
for (id in ensemble_ids) {
    gene_name = ensemble_to_gene_name(id)
    print(paste(id, gene_name))
}