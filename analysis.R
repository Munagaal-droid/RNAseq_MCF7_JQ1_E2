#!/usr/bin/env Rscript

user_library <- Sys.getenv(
  "R_LIBS_USER",
  unset = file.path(Sys.getenv("HOME"), ".local", "R", "library")
)
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
})

root_dir <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root_dir, "data")
results_dir <- file.path(root_dir, "results")
plot_dir <- file.path(results_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(data_dir, "gene_counts.csv")
metadata_path <- file.path(data_dir, "sample_metadata.csv")

if (!file.exists(counts_path)) {
  stop(
    "Missing data/gene_counts.csv. Run Rscript scripts/prepare_counts.R first.",
    call. = FALSE
  )
}

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
count_data <- read.csv(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
gene_ids <- count_data$gene_id
counts <- as.matrix(count_data[, setdiff(names(count_data), "gene_id"), drop = FALSE])
rownames(counts) <- gene_ids
storage.mode(counts) <- "numeric"

if (anyDuplicated(rownames(counts))) {
  stop("Gene identifiers in data/gene_counts.csv must be unique.", call. = FALSE)
}
if (!setequal(colnames(counts), metadata$sample_id)) {
  stop("Count-matrix columns do not match sample_metadata.csv.", call. = FALSE)
}
counts <- counts[, metadata$sample_id, drop = FALSE]
if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
  stop("The gene-count matrix contains missing or invalid values.", call. = FALSE)
}

metadata$treatment <- factor(metadata$treatment, levels = c("Control", "JQ1"))
metadata$estradiol <- factor(metadata$estradiol, levels = c("Veh", "E2"))
design <- model.matrix(~ estradiol + treatment, data = metadata)
if (qr(design)$rank != ncol(design)) {
  stop("The treatment design matrix is not full rank.", call. = FALSE)
}

write.csv(metadata, file.path(results_dir, "sample_metadata.csv"), row.names = FALSE)
write.csv(
  data.frame(
    gene_id = rownames(counts),
    mean_scaled_count = rowMeans(counts),
    median_scaled_count = apply(counts, 1, median),
    stringsAsFactors = FALSE
  ),
  file.path(results_dir, "count_summary.csv"),
  row.names = FALSE
)

dge <- DGEList(counts = counts, genes = data.frame(gene_id = rownames(counts)))
counts_per_million <- cpm(dge)
keep <- rowSums(counts_per_million > 1) >= 2
if (!any(keep)) {
  stop("No genes passed the CPM > 1 in at least two samples filter.", call. = FALSE)
}
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)

png(file.path(plot_dir, "scaled_count_boxplots.png"), width = 1400, height = 900, res = 140)
boxplot(log2(counts + 1), las = 2, ylab = "Log2 scaled counts + 1",
        main = "Recount3 gene counts by sample")
dev.off()

png(file.path(plot_dir, "library_sizes.png"), width = 1400, height = 900, res = 140)
barplot(dge$samples$lib.size, names.arg = colnames(dge), las = 2,
        ylab = "Estimated library size", main = "Library sizes")
dev.off()

png(file.path(plot_dir, "log_cpm_boxplots.png"), width = 1400, height = 900, res = 140)
boxplot(log_cpm, las = 2, ylab = "Log2 counts per million",
        main = "TMM-normalized log-CPM")
abline(h = median(log_cpm), col = "blue")
dev.off()

png(file.path(plot_dir, "mds.png"), width = 1200, height = 900, res = 140)
plotMDS(dge, labels = colnames(dge), col = as.integer(metadata$treatment),
        main = "MDS of TMM-normalized samples")
legend("topright", legend = levels(metadata$treatment),
       col = seq_along(levels(metadata$treatment)), pch = 16)
dev.off()

gene_variances <- apply(log_cpm, 1, var)
top_genes <- names(sort(gene_variances, decreasing = TRUE))[
  seq_len(min(500L, length(gene_variances)))
]
write.csv(
  data.frame(gene_id = top_genes, variance = gene_variances[top_genes]),
  file.path(results_dir, "top_variable_genes.csv"),
  row.names = FALSE
)

png(file.path(plot_dir, "top_variable_genes_heatmap.png"),
    width = 1400, height = 1400, res = 140)
heatmap(log_cpm[top_genes, , drop = FALSE],
        main = sprintf("Top %d most variable genes", length(top_genes)),
        scale = "row")
dev.off()

pca <- prcomp(t(log_cpm), scale. = TRUE)
pca_data <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  treatment = metadata$treatment,
  estradiol = metadata$estradiol
)
write.csv(pca_data, file.path(results_dir, "pca_scores.csv"), row.names = FALSE)
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = treatment, shape = estradiol,
                                 label = sample_id)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  theme_bw() +
  labs(title = "PCA of TMM-normalized log-CPM")
ggsave(file.path(plot_dir, "pca.png"), pca_plot, width = 9, height = 7, dpi = 150)

dge <- estimateDisp(dge, design, robust = TRUE)
png(file.path(plot_dir, "biological_coefficient_of_variation.png"),
    width = 1200, height = 900, res = 140)
plotBCV(dge, main = "Biological coefficient of variation")
dev.off()

fit <- glmFit(dge, design)
treatment_coefficient <- match("treatmentJQ1", colnames(design))
if (is.na(treatment_coefficient)) {
  stop("Could not identify the JQ1 coefficient in the design matrix.", call. = FALSE)
}
lrt <- glmLRT(fit, coef = treatment_coefficient)
edgeR_result <- topTags(
  lrt, n = Inf, adjust.method = "BH", sort.by = "PValue"
)$table
edgeR_result$gene_id <- rownames(edgeR_result)

ensembl_ids <- sub("\\.[0-9]+$", "", edgeR_result$gene_id)
gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(ensembl_ids),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
entrez_ids <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(ensembl_ids),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
edgeR_result$gene_name <- unname(gene_symbols[ensembl_ids])
edgeR_result$entrez_id <- unname(entrez_ids[ensembl_ids])
edgeR_result <- edgeR_result[
  , c("gene_id", "gene_name", "entrez_id",
      setdiff(names(edgeR_result), c("gene_id", "gene_name", "entrez_id"))),
  drop = FALSE
]
write.csv(edgeR_result, file.path(results_dir, "edgeR_results.csv"), row.names = FALSE)

degs <- edgeR_result[!is.na(edgeR_result$FDR) & edgeR_result$FDR < 0.05, , drop = FALSE]
upregulated <- degs[degs$logFC > 0, , drop = FALSE]
downregulated <- degs[degs$logFC < 0, , drop = FALSE]
write.csv(upregulated, file.path(results_dir, "upregulated_genes.csv"), row.names = FALSE)
write.csv(downregulated, file.path(results_dir, "downregulated_genes.csv"), row.names = FALSE)

de_gene_ids <- rownames(lrt)[as.logical(decideTestsDGE(lrt, p = 0.05))]
png(file.path(plot_dir, "smear_plot.png"), width = 1200, height = 900, res = 140)
plotSmear(lrt, de.tags = de_gene_ids, main = "JQ1 versus control")
abline(h = c(-1, 1), col = "blue")
dev.off()

if (nrow(degs) > 0L) {
  tested_entrez <- unique(na.omit(edgeR_result$entrez_id))
  deg_entrez <- unique(na.omit(degs$entrez_id))
  if (length(deg_entrez) == 0L) {
    stop("No significant genes mapped to human Entrez IDs for KEGG.", call. = FALSE)
  }
  kegg_links <- KEGGREST::keggLink("pathway", "hsa")
  kegg_pathways <- KEGGREST::keggList("pathway", "hsa")
  gene_pathways <- data.frame(
    entrez_id = sub("^hsa:", "", names(kegg_links)),
    pathway_id = sub("^path:", "", unname(kegg_links)),
    stringsAsFactors = FALSE
  )
  tested_entrez <- intersect(tested_entrez, gene_pathways$entrez_id)
  deg_entrez <- intersect(deg_entrez, tested_entrez)
  if (length(deg_entrez) == 0L) {
    stop("Significant genes did not map to KEGG pathways.", call. = FALSE)
  }
  pathway_genes <- split(gene_pathways$entrez_id, gene_pathways$pathway_id)
  tested_gene_count <- length(tested_entrez)
  deg_gene_count <- length(deg_entrez)
  pathway_rows <- lapply(names(pathway_genes), function(pathway_id) {
    members <- intersect(unique(pathway_genes[[pathway_id]]), tested_entrez)
    if (length(members) == 0L) {
      return(NULL)
    }
    overlap <- intersect(deg_entrez, members)
    data.frame(
      ID = pathway_id,
      Description = unname(kegg_pathways[paste0("path:", pathway_id)]),
      GeneRatio = sprintf("%d/%d", length(overlap), deg_gene_count),
      BgRatio = sprintf("%d/%d", length(members), tested_gene_count),
      pvalue = phyper(
        length(overlap) - 1L,
        length(members),
        tested_gene_count - length(members),
        deg_gene_count,
        lower.tail = FALSE
      ),
      geneID = paste(overlap, collapse = "/"),
      Count = length(overlap),
      stringsAsFactors = FALSE
    )
  })
  pathway_results <- do.call(rbind, pathway_rows)
  if (is.null(pathway_results) || nrow(pathway_results) == 0L) {
    stop("No KEGG pathways overlapped the significant genes.", call. = FALSE)
  }
  pathway_results$p.adjust <- p.adjust(pathway_results$pvalue, method = "BH")
  pathway_results <- pathway_results[order(pathway_results$p.adjust), ]
  significant_pathways <- pathway_results[
    !is.na(pathway_results$p.adjust) & pathway_results$p.adjust < 0.05, ,
    drop = FALSE
  ]
  write.csv(pathway_results, file.path(results_dir, "kegg_enrichment.csv"),
            row.names = FALSE)
  write.csv(significant_pathways, file.path(results_dir, "significant_pathways.csv"),
            row.names = FALSE)
  if (nrow(significant_pathways) > 0L) {
    plot_data <- head(significant_pathways, 10L)
    pathway_plot <- ggplot(
      plot_data,
      aes(x = Count, y = reorder(Description, Count), color = -log10(p.adjust),
          size = Count)
    ) +
      geom_point() +
      theme_bw() +
      labs(x = "Gene count", y = NULL, color = "-log10 adjusted p-value",
           title = "KEGG pathway enrichment")
    ggsave(file.path(plot_dir, "kegg_pathways.png"), pathway_plot,
           width = 10, height = 7, dpi = 150)
  }
} else {
  write.csv(data.frame(), file.path(results_dir, "kegg_enrichment.csv"),
            row.names = FALSE)
  write.csv(data.frame(), file.path(results_dir, "significant_pathways.csv"),
            row.names = FALSE)
  message("No genes met the FDR < 0.05 threshold; KEGG enrichment was skipped.")
}

session_info <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
writeLines(session_info, file.path(results_dir, "sessionInfo.txt"))
message(
  sprintf(
    "Analysis complete: %d genes tested; %d DEGs at FDR < 0.05.",
    nrow(edgeR_result), nrow(degs)
  )
)
