source("../deg/parse_counts.R")

genetype_lookup_filename = "../artifacts/genetype_lookup.txt"
genetype_lookup_df = get_genetype_lookup(genetype_lookup_filename)


ensemble_to_gene_name <- function(ensembl_id) {
  gene_name <- genetype_lookup_df$gene_name[genetype_lookup_df$gene_id == ensembl_id]
  if (length(gene_name) == 1) {
    return(gene_name[1])
  }
  else if (length(gene_name) > 1) {
    warning(paste("Multiple gene names found for Ensembl ID:", ensembl_id, "Returning the first one:", gene_name[1]))
    return(gene_name[1])
  } 
  else {
    stop(paste("No gene name found for Ensembl ID:", ensembl_id))
  }
}

gene_name_to_ensemble <- function(gene_name) {
    ensembl_id <- genetype_lookup_df$gene_id[genetype_lookup_df$gene_name == gene_name]
    if (length(ensembl_id) == 1) {
        return(ensembl_id[1])
    } 
    else if (length(ensembl_id) > 1) {
        warning(paste("Multiple Ensembl IDs found for gene name:", gene_name, "Returning the first one:", ensembl_id[1]))
        return(ensembl_id[1])
    }
    else {
        stop(paste("No Ensembl ID found for gene name:", gene_name))
    }
}
