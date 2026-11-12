#!/usr/bin/env Rscript

library_dir <- Sys.getenv(
  "R_LIBS_USER",
  unset = file.path(Sys.getenv("HOME"), ".local", "R", "library")
)
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = library_dir)
}

cran_packages <- c("ggplot2")
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org", lib = library_dir)
}

bioc_packages <- c("edgeR", "AnnotationDbi", "org.Hs.eg.db", "KEGGREST")
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0L) {
  BiocManager::install(missing_bioc, lib = library_dir, ask = FALSE, update = FALSE)
}
