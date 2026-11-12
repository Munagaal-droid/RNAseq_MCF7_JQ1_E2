#!/usr/bin/env Rscript

root_dir <- normalizePath(".", mustWork = TRUE)
data_dir <- file.path(root_dir, "data")
raw_dir <- file.path(data_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- paste0(
  "https://recount-opendata.s3.amazonaws.com/recount3/release/human/",
  "data_sources/sra/"
)
count_url <- paste0(
  base_url,
  "gene_sums/36/SRP040136/sra.gene_sums.SRP040136.G029.gz"
)
qc_url <- paste0(
  base_url,
  "metadata/36/SRP040136/sra.recount_qc.SRP040136.MD.gz"
)
count_path <- file.path(raw_dir, "SRP040136_G029_gene_sums.tsv.gz")
qc_path <- file.path(raw_dir, "SRP040136_recount_qc.tsv.gz")

download_if_missing <- function(url, destination) {
  if (file.exists(destination) && file.info(destination)$size > 0) {
    return(invisible(destination))
  }
  partial_path <- paste0(destination, ".part")
  status <- utils::download.file(
    url, partial_path, mode = "wb", method = "libcurl", quiet = FALSE
  )
  if (status != 0L || !file.rename(partial_path, destination)) {
    stop("Could not download data from ", url, call. = FALSE)
  }
  invisible(destination)
}

download_if_missing(count_url, count_path)
download_if_missing(qc_url, qc_path)

metadata <- read.csv(
  file.path(data_dir, "sample_metadata.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
coverage_table <- read.delim(
  gzfile(count_path),
  skip = 2,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
qc <- read.delim(
  gzfile(qc_path),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(metadata$run_accession %in% names(coverage_table))) {
  stop("The recount3 matrix is missing one or more requested SRA runs.", call. = FALSE)
}
qc_rows <- match(metadata$run_accession, qc$external_id)
auc <- qc$bc_auc.all_reads_all_bases[qc_rows]
if (anyNA(auc) || any(!is.finite(auc)) || any(auc <= 0)) {
  stop("Missing or invalid recount3 AUC values for one or more samples.", call. = FALSE)
}

coverage <- as.matrix(coverage_table[, metadata$run_accession, drop = FALSE])
rownames(coverage) <- coverage_table$gene_id
storage.mode(coverage) <- "numeric"
if (anyNA(coverage) || any(!is.finite(coverage)) || any(coverage < 0)) {
  stop("The recount3 coverage matrix contains invalid values.", call. = FALSE)
}

target_size <- 4e7
scale_factors <- target_size / auc
scaled_counts <- round(sweep(coverage, 2, scale_factors, `*`))
colnames(scaled_counts) <- metadata$sample_id
count_output <- data.frame(
  gene_id = rownames(scaled_counts),
  scaled_counts,
  check.names = FALSE
)
write.csv(
  count_output,
  file.path(data_dir, "gene_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  data.frame(
    sample_id = metadata$sample_id,
    run_accession = metadata$run_accession,
    auc = auc,
    scale_factor = scale_factors
  ),
  file.path(data_dir, "library_scaling_factors.csv"),
  row.names = FALSE
)

message(
  sprintf(
    "Prepared %d genes across %d samples using recount3 AUC scaling.",
    nrow(scaled_counts), ncol(scaled_counts)
  )
)
