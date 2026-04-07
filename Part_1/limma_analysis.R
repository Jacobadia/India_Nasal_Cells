library(edgeR)
library(limma)
library(readr)
library(tidyverse)
library(DESeq2)
library(sva)
library(org.Hs.eg.db)

setwd("/India_Nasal_Cells/data")

# Read and format the featureCounts table properly
counts <- read_tsv("gene_counts.tsv")
column_to_rownames("Geneid")

setwd("India_Nasal_Cells/Part_1")

# counts <- counts %>%
#   dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
#   column_to_rownames("Geneid") # Moves Geneid to row names so the matrix is strictly numeric

# Extract groups and create the DGEList
cols <- colnames(counts)
group <- factor(substr(cols, 12, 13))

d0 <- DGEList(counts = counts, group = group)

metadata <- read_tsv("metadata.tsv")
sex <- metadata %>%
  dplyr::select(Sex) %>%
  pull() %>%
  factor()

# Filter lowly expressed genes
keep <- filterByExpr(d0, group = group)
d <- d0[keep, , keep.lib.sizes = FALSE]

# Normalize
d <- calcNormFactors(d)

plotMDS(d, col = as.numeric(d$samples$group))

# Voom
mm <- model.matrix(~0 + group)

# Heatmap of model matrix (optional)
# library(pheatmap)
# pheatmap(mm, cluster_rows = FALSE, cluster_cols = FALSE)

voom.y.d <- voom(d, mm, plot = TRUE)

# Fit linear model and Contrasts
fit <- lmFit(voom.y.d, mm)

contr <- makeContrasts(LA_vs_LB = groupLA - groupLB, levels = colnames(mm))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

# Extract Results
top.table <- topTable(tmp, sort.by = "P", n = Inf)

# View DEGs
DEGs <- top.table %>% 
  arrange(adj.P.Val) #%>% 
  #filter(adj.P.Val < 0.5)

print(paste("Significant DEGs found:", nrow(DEGs)))
head(DEGs, 10)

clean_ids <- gsub("\\..*", "", rownames(DEGs))

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = clean_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

DEGs$Gene_Symbol <- gene_symbols

DEGs_final <- DEGs %>%
  rownames_to_column(var = "Ensembl_ID") %>%
  dplyr::select(Gene_Symbol, Ensembl_ID, everything())

head(DEGs_final)

#initial LIMMA Analysis concludes here


# ### UGANDA: Fisher Method for external dataset validation ###
# setwd("/home/bk249/groups/grp_tb/processing_scripts/DATA/UGANDA")

# #### Preparing data  ####
# counts_data <- read.csv("nasalcounts.csv", header=TRUE, row.names = 1)
# coldata <- read.csv("nasalcoldata.csv", header=TRUE, row.names = 1)
# coldata <- coldata %>%
#   dplyr::mutate(
#     status = factor(status, levels = c("control", "case")),
#     sex = factor(sex, levels = c("male", "female"))
#   )
# str(coldata)

# counts_data <- as.matrix(counts_data)
# all(rownames(coldata) %in% colnames(counts_data))
# all(rownames(coldata) == colnames(counts_data))
# all(colnames(counts_data) %in% rownames(coldata))
# all(colnames(counts_data) == rownames(coldata))


# #sex correction
# batch <- factor(coldata$sex)
# sample_group <- factor(coldata$status)
# counts_corrected <- ComBat_seq(as.matrix(counts_data),
#                                batch=batch,
#                                group=sample_group)

# #### DESeq2  ####
# dds <- DESeqDataSetFromMatrix(countData = counts_corrected,
#                               colData = coldata,
#                               design = ~ status)
# dds
# mcols(dds)
# dds$status <- relevel(dds$status, ref = "control")
# levels(dds$status)
# dds <- DESeq(dds)
# res05 <- results(dds, alpha=0.05)
# resOrdered_05 <- res05[order(res05$pvalue),]
# de <- as.data.frame(resOrdered_05)
# signif_de <- subset(de, padj < 0.05 & abs(log2FoldChange) >1.0)
# sig_genes <- row.names(signif_de)
# sig_genes
# sig_genes_nasal <- sig_genes
# table_counts_normalized <- counts(dds, normalized=TRUE)


# #### Parse files ####
# # normalised gene expression table
# em <- read.csv("data/em_nasal.csv", header=TRUE, row.names = 1)


# # differential expression table
# de <- read.csv("data/de_nasal.csv", header=TRUE, row.names = 1)
# res_nasal <- de
# de <- de[,-c(1,3,4)]
# colnames(de) <- c("log2fold", "p", "p.adj")

# # Create a fresh copy so 'df' remains untouched
# DEGs_mapped <- DEGs

# # Identify where your IDs are. 
# ids_to_map <- rownames(DEGs_mapped)

# # Clean the IDs/remove the .#######
# clean_ids <- gsub("\\..*", "", ids_to_map)


# # Perform the mapping
# library(org.Hs.eg.db)
# gene_symbols <- mapIds(org.Hs.eg.db,
#                        keys = clean_ids,
#                        column = "SYMBOL",
#                        keytype = "ENSEMBL",
#                        multiVals = "first")

# # Add the new symbols as a column in the NEW dataframe
# DEGs_mapped$Gene_Symbol <- gene_symbols

# # View the results
# head(DEGs_mapped)


# df_clean <- de %>% 
#   rownames_to_column(var = "Gene_Symbol")

# merged_data <- DEGs_mapped %>%
#   full_join(df_clean, by = "Gene_Symbol")

# merged_data_clean <- merged_data%>%
#   dplyr::select("Gene_Symbol", "P.Value", "adj.P.Val", "p","p.adj")

# merged_data_clean <- merged_data_clean %>%
#   dplyr::rename("Nis_P.adj" = "p.adj")

# merged_data_clean <- merged_data_clean %>%
#   dplyr::rename("Nis_raw.P" = "p")

# original_rows <- nrow(merged_data_clean)
# print(paste("Original number of rows:", original_rows))

# merged_data_clean <- merged_data_clean %>%
#   drop_na()

# #----
# # Fisher's method

# p1 <- merged_data_clean$P.Value
# p2 <- merged_data_clean$Nis_raw.P

# p1[p1 == 0] <- .Machine$double.xmin
# p2[p2 == 0] <- .Machine$double.xmin

# chi_sq_stats <- -2* (log(p1)+log(p2))

# merged_data_clean$Fisher_P <- pchisq(chi_sq_stats, df=4, lower.tail = FALSE)
# merged_data_clean$Fisher_Adj_P <- p.adjust(merged_data_clean$Fisher_P, method = "BH")
# merged_data_clean <- merged_data_clean %>%dplyr::select(Gene_Symbol,P.Value,Nis_raw.P,Fisher_P,Fisher_Adj_P)
# View(merged_data_clean)

# threshold = 0.05

# filtered_merged_clean <- merged_data_clean %>%
#   filter(merged_data_clean$Fisher_Adj_P<threshold)

# View(filtered_merged_clean)