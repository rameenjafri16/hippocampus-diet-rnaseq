# BINF6110 - Assignment 2
# Differential Expression Analysis: Dietary Restriction Effects on Mouse Hippocampal Transcriptome
# Dataset: GSE111778 - Wahl et al. (2018), Cell Reports
# Author: Rameen Jafri
#*********************************************************************

# LOADING LIBRARIES  ------
library(DESeq2)
library(tximport)
library(apeglm)
library(org.Mm.eg.db)    
library(clusterProfiler)
#library(enrichplot)    
library(pheatmap)
library(tidyverse)
library(Biostrings)
library(RColorBrewer)
library(ggrepel)
library(ggplotify)
#library(patchwork)
#library(knitr)
select <- AnnotationDbi::select

# PATH SETUP ------
setwd("~/BINF_DESKTOP/BINF6110/assignment02")

# LOADING METADATA ------
# Load sample information - one row per sample
metadata <- read.csv("data/metadata.csv", row.names = 1)

# Make diet a factor with Chow as the reference level
# This means all comparisons will be made against Chow
metadata$diet <- factor(metadata$diet, levels = c("Chow", "ChowCR", "P5", "P10", "P15"))
metadata$sex <- factor(metadata$sex)

# Check metadata looks right
print(metadata)
table(metadata$diet, metadata$sex)

# BUILDING TX2GENE MAPPING ------
# Salmon outputs transcript-level counts (ENSMUST IDs)
# DESeq2 works at the gene level (ENSMUSG IDs)
# tx2gene maps transcripts to genes so tximport can summarize them
# This is important because one gene can have multiple transcripts (alternative splicing)

# Read the FASTA file we used to build the Salmon index
fasta <- readDNAStringSet("reference/Mus_musculus.GRCm38.cdna.all.fa.gz")
headers <- names(fasta)

# Each header looks like:
# ENSMUST00000177564.1 cdna chromosome:GRCm38:... gene:ENSMUSG00000096176.1 ...
# We extract the transcript ID (before the first space) and gene ID (after "gene:")

# Extract transcript ID (everything before the first space)
tx_id <- sub(" .*", "", headers)

# Extract gene ID (the ENSMUSG... part after "gene:")
gene_id <- sub(".*gene:(ENSMUSG[0-9]+\\.[0-9]+).*", "\\1", headers)

# Build tx2gene dataframe with two columns: transcript ID and gene ID
tx2gene_ensembl <- data.frame(tx = tx_id, gene = gene_id)

head(tx2gene_ensembl)
nrow(tx2gene_ensembl)

# IMPORTING SALMON OUTPUT WITH TXIMPORT -----
# Build paths to all 30 quant.sf files
files <- file.path("salmon_output", rownames(metadata), "quant.sf")
names(files) <- rownames(metadata)

# Check all files exist before importing
all(file.exists(files))
length(files)

# Import transcript-level counts and summarize to gene level using tx2gene mapping
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene_ensembl,
                ignoreTxVersion = FALSE)  # Keep version numbers, tx2gene has them

# Check how many there are
nrow(txi$counts)
head(rownames(txi$counts))

# Do the tx IDs in tx2gene match the IDs in quant.sf?
head(rownames(txi$counts))

# CREATING DESEQ2 DATASET -----
# ~ sex + diet means: "account for sex differences, then look at diet effects"
dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ sex + diet)

# Sanity check: visualise normalised counts for Gpr17 (ENSMUSG00000052229.5)
# Top DEG across all protein restriction groups to confirm dose-dependent upregulation (coloured by sex) 
plotCounts(dds, gene="ENSMUSG00000052229.5", intgroup="diet")

# Get count data
gpr17_data <- plotCounts(
  dds,
  gene = "ENSMUSG00000052229.5",
  intgroup = c("diet","sex"),
  returnData = TRUE)

p_gpr17 <- ggplot(gpr17_data, aes(x = diet, y = count, color = sex)) +
  geom_point(size = 3, position = position_jitter(width = 0.1)) +
  scale_color_manual(values = c(male = "steelblue", female = "pink3")) +
  scale_x_discrete(limits = c("Chow","ChowCR","P5","P10","P15")) +
  labs(title = "Gpr17 — Normalized Counts by Diet Group",
       x = "Diet Group", y = "Normalized Count", color = "Sex") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p_gpr17
#ggsave("results/figures/Gpr17_counts.png",plot = p_gpr17,width = 7,height = 6,dpi = 300)

cat("Genes before filtering:", nrow(dds), "\n")

# Pre-filter: remove genes with very low counts across all samples
# Keeps genes with at least 10 reads in at least 6 samples (smallest group size)
keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep, ]
cat("Genes remaining after filtering:", nrow(dds), "\n")

# Run DESeq2 (performs normalization, dispersion estimation, and statistical testing)
dds <- DESeq(dds)

# Check what comparisons are available
resultsNames(dds)

# EXTRACTING RESULTS - EACH DIET VS CHOW -----
# lfcShrink improves fold change estimates for lowly expressed genes
# using apeglm shrinkage (Zhu et al. 2018)

# ChowCR vs Chow
res_CR <- lfcShrink(dds, coef = "diet_ChowCR_vs_Chow", type = "apeglm")

# 5% protein vs Chow
res_P5 <- lfcShrink(dds, coef = "diet_P5_vs_Chow", type = "apeglm")

# 10% protein vs Chow
res_P10 <- lfcShrink(dds, coef = "diet_P10_vs_Chow", type = "apeglm")

# 15% protein vs Chow
res_P15 <- lfcShrink(dds, coef = "diet_P15_vs_Chow", type = "apeglm")

# Summary of significant genes in each comparison
cat("\n--- ChowCR vs Chow ---\n"); summary(res_CR, alpha = 0.05)
cat("\n--- 5% Protein vs Chow ---\n"); summary(res_P5, alpha = 0.05)
cat("\n--- 10% Protein vs Chow ---\n"); summary(res_P10, alpha = 0.05)
cat("\n--- 15% Protein vs Chow ---\n"); summary(res_P15, alpha = 0.05)

# MA PLOTS -----
# MA plots show fold change vs mean expression
# Good for checking normalization worked correctly
# need to convert base plots to ggplot objects
p1 <- as.ggplot(~plotMA(res_CR,  main="ChowCR vs Chow",  ylim=c(-3,3), alpha=0.05))
p2 <- as.ggplot(~plotMA(res_P5,  main="5% Protein vs Chow", ylim=c(-3,3), alpha=0.05))
p3 <- as.ggplot(~plotMA(res_P10, main="10% Protein vs Chow", ylim=c(-3,3), alpha=0.05))
p4 <- as.ggplot(~plotMA(res_P15, main="15% Protein vs Chow", ylim=c(-3,3), alpha=0.05))

ma_combined <- (p1 | p2) / (p3 | p4) # Combine into 2x2 grid
ma_combined

#ggsave("results/figures/MA_plots.png", plot = ma_combined, width = 12, height = 10, dpi = 300)

# Top 10 differentially regulated genes when comparing each group to 19% protein (chow diet)
get_top10 <- function(res, comparison_name) {
  
  res_df <- as.data.frame(res)
  res_df$gene_clean <- sub("\\..*", "", rownames(res_df))
  
  # Add gene symbols
  res_df$symbol <- mapIds(org.Mm.eg.db,
                          keys = res_df$gene_clean,
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
  
  # Get top 10 by adjusted p-value - avoid select() conflict
  top10 <- res_df[!is.na(res_df$padj), ]
  top10 <- top10[order(top10$padj), ]
  top10 <- head(top10, 10)
  top10 <- top10[, c("symbol", "log2FoldChange", "padj")]
  
  cat("\n--- Top 10 DEGs:", comparison_name, "---\n")
  print(top10)
  return(top10)}

top10_CR  <- get_top10(res_CR,  "CR vs Chow")
top10_P5  <- get_top10(res_P5,  "5% Protein vs Chow")
top10_P10 <- get_top10(res_P10, "10% Protein vs Chow")
top10_P15 <- get_top10(res_P15, "15% Protein vs Chow")

# Printing in wide format to match the paper's Table S2 style ------
top10_wide <- data.frame(
  Rank = 1:10,
  CR_Gene     = top10_CR$symbol,
  CR_log2FC   = round(top10_CR$log2FoldChange, 3),
  CR_padj     = formatC(top10_CR$padj, format="e", digits=2),
  P5_Gene     = top10_P5$symbol,
  P5_log2FC   = round(top10_P5$log2FoldChange, 3),
  P5_padj     = formatC(top10_P5$padj, format="e", digits=2),
  P10_Gene    = top10_P10$symbol,
  P10_log2FC  = round(top10_P10$log2FoldChange, 3),
  P10_padj    = formatC(top10_P10$padj, format="e", digits=2),
  P15_Gene    = top10_P15$symbol,
  P15_log2FC  = round(top10_P15$log2FoldChange, 3),
  P15_padj    = formatC(top10_P15$padj, format="e", digits=2))

print(top10_wide)
#write.csv(top10_wide, "results/top10_DEGs_wide.csv", row.names=FALSE)

#DEG count summary ------
deg_summary <- data.frame(
  Comparison = c("CR vs Chow", "5% Protein vs Chow", 
                 "10% Protein vs Chow", "15% Protein vs Chow"),
  Total_DEGs = c(
    sum(res_CR$padj < 0.05, na.rm=TRUE),
    sum(res_P5$padj < 0.05, na.rm=TRUE),
    sum(res_P10$padj < 0.05, na.rm=TRUE),
    sum(res_P15$padj < 0.05, na.rm=TRUE)
  ),
  Upregulated = c(
    sum(res_CR$padj < 0.05 & res_CR$log2FoldChange > 0, na.rm=TRUE),
    sum(res_P5$padj < 0.05 & res_P5$log2FoldChange > 0, na.rm=TRUE),
    sum(res_P10$padj < 0.05 & res_P10$log2FoldChange > 0, na.rm=TRUE),
    sum(res_P15$padj < 0.05 & res_P15$log2FoldChange > 0, na.rm=TRUE)
  ),
  Downregulated = c(
    sum(res_CR$padj < 0.05 & res_CR$log2FoldChange < 0, na.rm=TRUE),
    sum(res_P5$padj < 0.05 & res_P5$log2FoldChange < 0, na.rm=TRUE),
    sum(res_P10$padj < 0.05 & res_P10$log2FoldChange < 0, na.rm=TRUE),
    sum(res_P15$padj < 0.05 & res_P15$log2FoldChange < 0, na.rm=TRUE)))

print(deg_summary)
#write.csv(deg_summary, "results/DEG_summary_table1.csv", row.names=FALSE)

# VARIANCE STABILIZING TRANSFORMATION FOR VISUALIZATION -----
# VST puts counts on a log-like scale suitable for visualization
# blind = FALSE means it uses the design to stabilize variance
vsd <- vst(dds, blind = FALSE)

# Compute sample-to-sample distances from VST counts -----
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Build annotation dataframe with diet and sex
annotation_df <- data.frame(
  Diet = metadata$diet,
  Sex = metadata$sex,
  row.names = rownames(metadata))

# Color scheme
ann_colors <- list(
  Diet = c(Chow = "gray50", ChowCR = "steelblue",
           P5 = "darkred", P10 = "tomato", P15 = "lightsalmon"),
  Sex = c(male = "skyblue", female = "pink"))

# Use a clean blue-white color scheme for distances
# (low distance = similar = dark blue, high distance = different = white)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
p_dist <- pheatmap(sampleDistMatrix,
                   clustering_distance_rows = sampleDists,
                   clustering_distance_cols = sampleDists,
                   col = colors,
                   annotation_row = annotation_df,
                   annotation_col = annotation_df,
                   annotation_colors = ann_colors,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   fontsize_row = 7,
                   fontsize_col = 7,
                   main = "Sample-to-Sample Distance Matrix")
#ggsave("results/figures/sample_distance_heatmap.png",plot = p_dist$gtable,width = 10,height = 9.5,dpi = 300)


# PCA PLOT -----
# PCA shows overall sample similarity - do samples cluster by diet group?
pca_data <- plotPCA(vsd, intgroup = c("diet", "sex"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = diet, shape = sex)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot: Hippocampal Gene Expression by Diet and Sex") +
  scale_color_manual(values = c("Chow"   = "gray50",
                                "ChowCR" = "steelblue",
                                "P5"     = "darkred",
                                "P10"    = "tomato",
                                "P15"    = "lightsalmon")) +
  theme_bw() +
  theme(legend.position = "right")

print(pca_plot)
#ggsave("results/figures/PCA_plot.png", pca_plot, width = 8, height = 6, dpi = 300)

# Look at the actual PCA coordinates to identify the outlier
pca_data <- plotPCA(vsd, intgroup=c("diet","sex"), returnData=TRUE)
print(pca_data)

# identify the outlier sample
pca_data[pca_data$PC2 > 10, ]

# HEATMAP OF TOP VARIABLE GENES -----
# Select top 50 most variable genes across all samples
top_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[top_var_genes, ]

# Scale each gene so expression levels are comparable across genes
mat <- t(scale(t(mat)))

# Annotation for columns (samples)
annotation_df <- metadata[, c("diet", "sex")]
colnames(annotation_df) <- c("Diet", "Sex")

# Color scheme for annotations
ann_colors <- list(
  Diet = c(Chow = "gray50", ChowCR = "steelblue",
           P5 = "darkred", P10 = "tomato", P15 = "lightsalmon"),
  Sex = c(male = "skyblue", female = "pink"))
p <- pheatmap(mat,
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              annotation_col = annotation_df,
              annotation_colors = ann_colors,
              show_rownames = FALSE,
              show_colnames = FALSE,
              main = "Top 50 Most Variable Genes")
#ggsave("results/figures/heatmap_top_50_variable_genes.png",plot = p$gtable,width = 10,height = 8)


# HEATMAP OF TOP DIFFERENTIALLY EXPRESSED GENES (CR vs Chow) -----
# Get top 50 significant DEGs from CR comparison, ordered by adjusted p-value
res_CR_df <- as.data.frame(res_CR) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50)

top_de_genes <- res_CR_df$gene

# Extract VST values for these genes
mat_de <- assay(vsd)[top_de_genes, ]

# Scale across rows
mat_de <- t(scale(t(mat_de)))

# Order columns by diet group so heatmap is easier to read
col_order <- order(factor(metadata$diet, 
                          levels = c("ChowCR", "P5", "P10", "P15", "Chow")))
mat_de <- mat_de[, col_order]
annotation_df_ordered <- annotation_df[col_order, ]

# Plot
p_cr <- pheatmap(mat_de,
                 cluster_rows = TRUE,
                 cluster_cols = FALSE,   # keep diet groups together
                 annotation_col = annotation_df_ordered,
                 annotation_colors = ann_colors,
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 main = "Top 50 DEGs: Caloric Restriction vs Chow")

#ggsave("results/figures/heatmap_top_CR_DEGs.png",plot = p_cr$gtable,width = 10,height = 8,dpi = 300)


# VOLCANO PLOTS - ONE PER COMPARISON -----
make_volcano <- function(res, title) {
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- na.omit(res_df)
  res_df$gene_clean <- sub("\\..*", "", res_df$gene)
  
  res_df$significant <- ifelse(res_df$padj < 0.05,
                               ifelse(res_df$log2FoldChange > 0, "Up", "Down"),
                               "Not Sig")
  
  res_df$symbol <- mapIds(org.Mm.eg.db,
                          keys = res_df$gene_clean,
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
  
  top_genes <- res_df %>%
    filter(significant != "Not Sig") %>%
    arrange(padj) %>%
    head(15)
  
  p <- ggplot(res_df, aes(x = log2FoldChange,
                          y = -log10(pvalue),
                          color = significant)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Down" = "steelblue",
                                  "Not Sig" = "gray70",
                                  "Up" = "darkred")) +
    geom_label_repel(data = top_genes,
                     aes(label = symbol),
                     size = 3,
                     max.overlaps = 20,
                     box.padding = 0.5,
                     show.legend = FALSE) +
    labs(x = "Log2 Fold Change",
         y = "-Log10 p-value",
         title = title,
         color = "Expression") +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed",
               color = "black",
               alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed",
               color = "black",
               alpha = 0.5) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)}
p_volcano_CR <- make_volcano(res_CR, "Caloric Restriction vs Chow")
p_volcano_CR
#ggsave("results/figures/volcano_CR.png", plot = p_volcano_CR, width = 9, height = 7, dpi = 300)

# FUNCTIONAL ENRICHMENT ANALYSIS - GSEA -----
# GSEA uses the full ranked gene list - detects pathway shifts without a significance cutoff
# Get unshrunken results just for ranking (keep shrunken for everything else)
res_CR_unshrunken  <- results(dds, contrast=c("diet","ChowCR","Chow"))
res_P5_unshrunken  <- results(dds, contrast=c("diet","P5","Chow"))
res_P10_unshrunken <- results(dds, contrast=c("diet","P10","Chow"))
res_P15_unshrunken <- results(dds, contrast=c("diet","P15","Chow"))

colnames(as.data.frame(res_CR_unshrunken))
run_GSEA <- function(res, comparison_name) {
  
  # Convert to dataframe explicitly and extract stat column
  res_df <- as.data.frame(res)
  res_df$gene_clean <- sub("\\..*", "", rownames(res_df))
  
  # Remove NAs from the stat column explicitly
  res_df <- res_df[!is.na(res_df$stat), ]
  
  cat(comparison_name, ": starting with", nrow(res_df), "genes\n")
  
  # Convert to Entrez IDs
  id_map <- bitr(res_df$gene_clean,
                 fromType = "ENSEMBL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
  
  # Merge by Ensembl ID
  res_df <- merge(res_df, id_map, by.x = "gene_clean", by.y = "ENSEMBL")
  
  # Create ranked named vector
  gene_list <- res_df$stat
  names(gene_list) <- res_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]
  
  cat(comparison_name, ": ranking", length(gene_list), "genes\n")
  
  gsea_result <- gseGO(geneList = gene_list,
                       OrgDb = org.Mm.eg.db,
                       ont = "BP",
                       minGSSize = 15,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       verbose = FALSE)
  
  cat(comparison_name, ": found", nrow(gsea_result), "enriched terms\n")
  return(gsea_result)}

gsea_CR  <- run_GSEA(res_CR_unshrunken,  "CR vs Chow")
gsea_P5  <- run_GSEA(res_P5_unshrunken,  "P5 vs Chow")
gsea_P10 <- run_GSEA(res_P10_unshrunken, "P10 vs Chow")
gsea_P15 <- run_GSEA(res_P15_unshrunken, "P15 vs Chow")

for (gsea_obj in list(list(gsea_CR,  "GSEA_CR",  "CR vs Chow"),
                      list(gsea_P5,  "GSEA_P5",  "5% Protein vs Chow"),
                      list(gsea_P10, "GSEA_P10", "10% Protein vs Chow"),
                      list(gsea_P15, "GSEA_P15", "15% Protein vs Chow"))) {
  
  if (!is.null(gsea_obj[[1]]) && nrow(gsea_obj[[1]]) > 0) {
    png(paste0("results/figures/", gsea_obj[[2]], ".png"), width=1100, height=900)
    print(dotplot(gsea_obj[[1]], showCategory=15, split=".sign") +
            facet_grid(.~.sign) +
            ggtitle(paste("GSEA GO BP:", gsea_obj[[3]])) +
            theme(plot.title = element_text(hjust=0.5, face="bold")))
    dev.off()
    system(paste0("open results/figures/", gsea_obj[[2]], ".png"))}}

cat("Analysis complete! Results saved to results\n")
