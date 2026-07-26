#!/usr/bin/env Rscript
# =============================================================================
# Install every R dependency used by this repository.
#
#   Rscript environment/install_R_packages.R
#
# Requires R >= 4.3 and network access. Bioconductor packages are installed
# via BiocManager; everything else comes from CRAN.
# =============================================================================

cran_packages <- c(
  "BiocManager",
  # tidyverse core (installed individually rather than as the metapackage so
  # a single failure does not block the rest)
  "dplyr", "ggplot2", "readr", "stringr", "tibble", "tidyr", "purrr",
  # plotting and layout
  "patchwork", "ggrepel", "pheatmap", "VennDiagram", "RColorBrewer",
  "scales", "futile.logger",
  # misc
  "ashr"
)

bioc_packages <- c(
  "DESeq2",
  "apeglm",
  "limma",
  "clusterProfiler",
  "AnnotationDbi",
  "org.Ce.eg.db",
  "celegans.db",
  "GO.db"
)

install_missing <- function(pkgs, installer, label) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) {
    message(sprintf("[ok] All %s packages already installed.", label))
    return(invisible(NULL))
  }
  message(sprintf("[install] %s: %s", label, paste(missing, collapse = ", ")))
  installer(missing)
}

repos <- getOption("repos")
if (is.null(repos[["CRAN"]]) || repos[["CRAN"]] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

install_missing(cran_packages,
                function(p) install.packages(p, dependencies = TRUE),
                "CRAN")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("BiocManager failed to install; cannot continue with Bioconductor packages.")
}

install_missing(bioc_packages,
                function(p) BiocManager::install(p, update = FALSE, ask = FALSE),
                "Bioconductor")

# ---- Report ------------------------------------------------------------------
all_pkgs <- c(setdiff(cran_packages, "BiocManager"), bioc_packages)
still_missing <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(still_missing)) {
  warning("Could not install: ", paste(still_missing, collapse = ", "))
} else {
  message("\n[done] All R dependencies present.")
}

cat("\n--- sessionInfo ---\n")
print(sessionInfo())
