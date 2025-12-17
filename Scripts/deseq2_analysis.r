############################################################
## Bulk RNA-seq: POS vs NEG (DESeq2 + GSEA)
## Author: Bhakti
## Description:
## Differential expression and pathway analysis with
## covariate adjustment (gender, age)
############################################################

## ----------------------------
## Libraries
## ----------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

## ----------------------------
## Create output directories
## ----------------------------
dir.create("results", showWarnings = FALSE)
dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("results/heatmaps", recursive = TRUE, showWarnings = FALSE)
dir.create("results/gsea", recursive = TRUE, showWarnings = FALSE)

## ----------------------------
## Load data
## ----------------------------
counts <- read.csv(here("data", "raw", "counts_filtered.csv"), row.names = 1)
meta   <- read.csv(here("data", "raw", "metadata.csv"), row.names = 1)


stopifnot(all(colnames(counts) == rownames(meta)))

## ----------------------------
## Metadata cleaning
## ----------------------------
meta$positivity <- tolower(trimws(meta$positivity))

meta$positivity[meta$positivity %in% c(
  "positive","pos","covid pos","covid+","ars-cov-2 positivity: pos"
)] <- "pos"

meta$positivity[meta$positivity %in% c(
  "negative","neg","covid neg","covid-","ars-cov-2 positivity: neg"
)] <- "neg"

meta$positivity <- factor(meta$positivity)
meta$gender     <- factor(meta$gender)
meta$batch      <- factor(meta$batch)

## ----------------------------
## Core DESeq2 + GSEA function
## ----------------------------
run_deseq2_gsea <- function(counts, meta, design, contrast,
                            prefix, gsea_ontology = "BP") {
  
  message("Running: ", prefix)
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = design
  )
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = contrast)
  res <- res[!is.na(res$padj), ]
  res <- res[order(res$padj), ]
  
  write.csv(res, file.path("results/deseq2",
                           paste0(prefix, "_DEG.csv")))
  
  ## Volcano
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  p_volcano <- ggplot(res_df,
                      aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(color = padj < 0.05), alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal()
  
  ggsave(file.path("results/deseq2",
                   paste0(prefix, "_volcano.png")),
         p_volcano, width = 6, height = 5)
  
  ## Heatmap (top 30 DEGs)
  sig <- res_df[res_df$padj < 0.05 &
                  abs(res_df$log2FoldChange) >= 1, ]
  
  if (nrow(sig) >= 2) {
    vsd <- vst(dds, blind = FALSE)
    top_genes <- head(sig$gene, 30)
    mat <- assay(vsd)[top_genes, ]
    
    png(file.path("results/heatmaps",
                  paste0(prefix, "_heatmap.png")))
    pheatmap(mat, scale = "row")
    dev.off()
  }
  
  ## GSEA
  entrez <- bitr(res_df$gene,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)
  
  merged <- merge(res_df, entrez,
                  by.x = "gene", by.y = "SYMBOL")
  
  gene_list <- merged$log2FoldChange
  names(gene_list) <- merged$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = gsea_ontology,
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    eps = 0
  )
  
  write.csv(gsea@result,
            file.path("results/gsea",
                      paste0(prefix, "_GSEA.csv")),
            row.names = FALSE)
  
  return(list(res = res, gsea = gsea))
}

## ----------------------------
## Unadjusted analysis
## ----------------------------
res_unadj <- run_deseq2_gsea(
  counts,
  meta,
  design   = ~ positivity,
  contrast = c("positivity", "pos", "neg"),
  prefix   = "pos_vs_neg_unadjusted"
)

## ----------------------------
## Adjusted analysis (gender + age)
## ----------------------------
meta$age <- as.numeric(meta$age)
meta_adj <- meta[!is.na(meta$age), ]
counts_adj <- counts[, rownames(meta_adj)]

res_adj <- run_deseq2_gsea(
  counts_adj,
  meta_adj,
  design   = ~ gender + age + positivity,
  contrast = c("positivity", "pos", "neg"),
  prefix   = "pos_vs_neg_adj_gender_age"
)
