# Figure 1 -- Acute cobalt (II) chloride exposure induces transcriptional
# responses in wildtype and arginine kinase loss-of-function mutants.
#
#
# Panels:
#   A: PCA of VST-transformed counts, with inter-experiment batch effect
#      removed via limma::removeBatchEffect (visualization only)
#   B: Knockout validation -- log2FC of target argk transcript in its
#      respective mutant vs N2 (untreated baseline)
#   C: Volcano plots -- argk-2(-/-) and argk-4(-/-) vs N2, in each treatment
#   G: Three-set Venn diagram of treatment-response DEGs across genotypes
#
# Pipeline decisions 
#   * argk-1 (MAH205) excluded -- incomplete knockdown in raw counts.
#
#   * PCA-only batch correction: limma::removeBatchEffect() is applied to
#     the VST matrix using experiment as the batch variable, with the
#     genotype*treatment design preserved. 
#   * DEGs for Panel D (Venn): BH-adjusted p < 0.05 AND |log2FC| > 1.
#     Panel C volcano highlighting uses the same thresholds. Both match
#     the published Figure 1 caption.



# 0. Setup 

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(stringr); library(purrr); library(ggplot2)
  library(DESeq2)
  library(limma)         # removeBatchEffect (for PCA only)
  library(ggrepel)
  library(VennDiagram)   # apt-available alternative to ggVennDiagram
  library(grid)
  library(patchwork)
  library(scales)
})
# VennDiagram writes a noisy log via futile.logger; silence it
suppressPackageStartupMessages(library(futile.logger))
flog.threshold(ERROR)

set.seed(2025)

# ---- Config ----
COUNTS_FILE   <- "data/RNASEQ61125.csv"
METADATA_FILE <- "data/Sample_Metadata_Table.csv"
OUT_DIR       <- "outputs/figures/Figure1"
DEG_DIR       <- file.path(OUT_DIR, "DEGs")

PADJ_CUTOFF   <- 0.05
LFC_CUTOFF    <- 1
MIN_COUNT_SUM <- 10

EXCLUDE_STRAINS  <- c("MAH205")

#          kept so select(-any_of(ANNOT_COLS)) remains a no-op gracefully.
ANNOT_COLS       <- c("ENTREZID", "ENSEMBL", "WORMBASE")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DEG_DIR, showWarnings = FALSE, recursive = TRUE)


# Style 

COL_UP   <- "#cb181d"
COL_DOWN <- "#2171b5"
COL_NS   <- "grey70"
COL_NA   <- "grey85"

COL_UNTR <- c(low = "#f7fbff", high = "#2171b5")
COL_TR   <- c(low = "#fff5f0", high = "#cb181d")

GROUP_COLORS <- c(
  "N2_untreated"     = "#3CB39C",
  "N2_treated"       = "#97C04F",
  "RB2060_untreated" = "#6FAEDB",
  "RB2060_treated"   = "#3FC5D8",
  "RB2598_untreated" = "#F2A6BC",
  "RB2598_treated"   = "#E66FA8"
)

GROUP_LABEL_EXPRS <- c(
  "N2_untreated"     = "'reference untreated'",
  "N2_treated"       = "'reference treated'",
  "RB2060_untreated" = "italic('argk-2')^'-/-'*' untreated'",
  "RB2060_treated"   = "italic('argk-2')^'-/-'*' treated'",
  "RB2598_untreated" = "italic('argk-4')^'-/-'*' untreated'",
  "RB2598_treated"   = "italic('argk-4')^'-/-'*' treated'"
)

GROUP_LEVELS <- c("N2_untreated",     "N2_treated",
                  "RB2060_untreated", "RB2060_treated",
                  "RB2598_untreated", "RB2598_treated")

GENOTYPE_LABELS <- c(
  "N2"     = "N2",
  "RB2060" = "argk-2(-/-)",
  "RB2598" = "argk-4(-/-)"
)

# [FIX 4] gene_id values updated with CELE_ prefix to match new counts file.
KO_TARGETS <- tibble(
  strain     = c("RB2060",       "RB2598"),
  gene_label = c("argk-2",       "argk-4"),
  gene_id    = c("CELE_W10C8.5", "CELE_F46H5.3")
)


# 1. Load and align data 

# Use robust CSV loader from shared utils (handles CR-only line endings)
# Source the utils file from one of these standard locations:
for (utils_path in c("scripts/_utils.R", "_utils.R", "../scripts/_utils.R")) {
  if (file.exists(utils_path)) { source(utils_path); break }
}
if (!exists("read_csv_robust")) {
  stop("Could not source scripts/_utils.R. Run this script from the repo root.")
}

raw_counts <- read_csv_robust(COUNTS_FILE, check.names = FALSE,
                              stringsAsFactors = FALSE)

# [FIX 1] Guard updated: new file uses "gene" not "gene_symbol"; no GENENAME.
stopifnot("gene" %in% colnames(raw_counts))

# [FIX 1] Renamed gene column; [FIX 2] biotype filter block removed entirely.
raw_counts <- raw_counts %>%
  filter(!is.na(gene),
         !duplicated(gene))

# Coerce sample columns to numeric (CR line endings can make them character)
sample_cols <- setdiff(colnames(raw_counts), c("gene", ANNOT_COLS))
for (cc in sample_cols) {
  raw_counts[[cc]] <- as.numeric(as.character(raw_counts[[cc]]))
}

cat("Genes loaded: ", nrow(raw_counts), "\n", sep = "")
cat("(Note: biotype filter not applied -- no GENENAME column in new counts file.)\n")

# [FIX 1] column renamed from "gene_symbol" to "gene"
counts <- raw_counts %>%
  column_to_rownames("gene") %>%
  select(-any_of(ANNOT_COLS))

coldata <- read_csv_robust(METADATA_FILE)
id_col_fig1 <- intersect(c("SampleLabel", "library_name", "Sample", "sample"),
                         colnames(coldata))[1]
if (is.na(id_col_fig1)) id_col_fig1 <- colnames(coldata)[1]
rownames(coldata) <- trimws(coldata[[id_col_fig1]])
coldata <- coldata[!coldata$genotype %in% EXCLUDE_STRAINS, , drop = FALSE]

# Harmonize sample naming convention (rename GEO library names to A-labels if needed)
counts_mat_temp <- as.matrix(counts)
counts_mat_temp <- harmonize_sample_names(counts_mat_temp, coldata,
                                          target_col = id_col_fig1,
                                          alt_col = setdiff(c("SampleLabel", "GEO_library_name"),
                                                            id_col_fig1)[1])
counts <- as.data.frame(counts_mat_temp, check.names = FALSE)

shared_samples <- intersect(colnames(counts), rownames(coldata))
if (length(shared_samples) == 0L) {
  stop("No overlapping sample names between counts and metadata.")
}
counts  <- counts[, shared_samples]
coldata <- coldata[shared_samples, , drop = FALSE]

counts <- as.matrix(counts)
mode(counts) <- "integer"

coldata$genotype   <- relevel(factor(coldata$genotype),   ref = "N2")
coldata$treatment  <- relevel(factor(coldata$treatment),  ref = "untreated")
coldata$experiment <- factor(coldata$experiment)

stopifnot(all(colnames(counts) == rownames(coldata)))

cat("\nLoaded ", nrow(counts), " genes x ", ncol(counts), " samples.\n", sep = "")
cat("Sample breakdown:\n")
print(table(coldata$genotype, coldata$treatment))


# 2. DESeq2 models 

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ genotype + treatment + experiment +
    genotype:treatment +
    genotype:experiment +
    treatment:experiment
)
dds <- dds[rowSums(counts(dds)) >= MIN_COUNT_SUM, ]
dds <- DESeq(dds)

make_within_dds <- function(tr_label) {
  cd <- coldata[coldata$treatment == tr_label, , drop = FALSE]
  cd$genotype <- relevel(droplevels(cd$genotype), ref = "N2")
  cc <- counts[, rownames(cd)]
  d  <- DESeqDataSetFromMatrix(countData = cc, colData = cd, design = ~ genotype)
  d  <- d[rowSums(counts(d)) >= MIN_COUNT_SUM, ]
  DESeq(d)
}
dds_un <- make_within_dds("untreated")
dds_tr <- make_within_dds("treated")

dds_n2 <- {
  cd <- coldata[coldata$genotype == "N2", , drop = FALSE]
  cd$treatment <- relevel(droplevels(cd$treatment), ref = "untreated")
  cc <- counts[, rownames(cd)]
  d  <- DESeqDataSetFromMatrix(countData = cc, colData = cd, design = ~ treatment)
  d  <- d[rowSums(counts(d)) >= MIN_COUNT_SUM, ]
  DESeq(d)
}


#  3. Panel A - PCA 

vsd <- vst(dds, blind = FALSE)

# Batch correction (VISUALIZATION ONLY).
mat <- assay(vsd)
mm  <- model.matrix(~ genotype * treatment,
                    data = as.data.frame(colData(vsd)))

assay(vsd) <- limma::removeBatchEffect(
  mat,
  batch  = colData(vsd)$experiment,
  design = mm
)

pcaData <- plotPCA(
  vsd,
  intgroup = c("genotype", "treatment", "experiment"),
  returnData = TRUE
)

percentVar <- round(100 * attr(pcaData, "percentVar"))

# FLIP PC2 FOR FIGURE ORIENTATION
#
# PCA eigenvector signs are arbitrary: PC2 and -PC2 represent the exact same
# principal component. This flip changes only the visual orientation so that
# Panel A matches the orientation used in the original manuscript figure.
# Variance explained, sample-to-sample distances, clustering, and all
# downstream analyses are unchanged.
# 
pcaData$PC2 <- -pcaData$PC2


cat("\nANOVA on principal components (after batch correction):\n")
cat("  PC1:\n")
print(
  summary(
    aov(PC1 ~ genotype * treatment + experiment, pcaData)
  )[[1]]
)

cat("  PC2:\n")
print(
  summary(
    aov(PC2 ~ genotype * treatment + experiment, pcaData)
  )[[1]]
)


pcaData <- pcaData %>%
  mutate(
    group = factor(
      paste(genotype, treatment, sep = "_"),
      levels = GROUP_LEVELS
    )
  )


panelA <- ggplot(
  pcaData,
  aes(
    x = PC1,
    y = PC2,
    colour = group
  )
) +
  geom_point(
    size = 3.5,
    alpha = 0.95
  ) +
  scale_colour_manual(
    name   = NULL,
    values = GROUP_COLORS,
    labels = parse(text = GROUP_LABEL_EXPRS[GROUP_LEVELS]),
    breaks = GROUP_LEVELS
  ) +
  guides(
    colour = guide_legend(
      ncol = 2,
      byrow = TRUE
    )
  ) +
  labs(
    x = sprintf(
      "PC1: %d%% variance",
      percentVar[1]
    ),
    y = sprintf(
      "PC2: %d%% variance",
      percentVar[2]
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    
    legend.background = element_rect(
      fill = alpha("white", 0.7),
      colour = NA
    ),
    
    legend.key.size = unit(
      0.4,
      "cm"
    ),
    
    legend.text = element_text(
      size = 10
    )
  )


ggsave(
  file.path(
    OUT_DIR,
    "Figure1_PanelA_PCA.pdf"
  ),
  panelA,
  width = 7.5,
  height = 5.5
)

# 4. Panel B -- Knockout validation 
# Thesis-style layout: two transcript-specific subpanels. Each subpanel shows
# both mutant genotypes vs N2 at untreated baseline, so the cognate KO and the
# non-cognate control comparison are visible side-by-side.

ko_rows <- list()
row_i <- 1L

for (target_i in seq_len(nrow(KO_TARGETS))) {
  gene_id    <- KO_TARGETS$gene_id[target_i]
  gene_label <- KO_TARGETS$gene_label[target_i]

  for (strain in c("RB2060", "RB2598")) {
    res <- as.data.frame(
      results(
        dds_un,
        contrast = c("genotype", strain, "N2"),
        alpha = PADJ_CUTOFF
      )
    )

    if (!gene_id %in% rownames(res)) {
      warning(sprintf("Target %s not found in dds_un results for %s.", gene_id, strain))
      next
    }

    ko_rows[[row_i]] <- tibble(
      target_gene = gene_label,
      strain      = strain,
      log2FC      = res[gene_id, "log2FoldChange"],
      padj        = res[gene_id, "padj"]
    )
    row_i <- row_i + 1L
  }
}

ko_lfc <- bind_rows(ko_rows)
if (nrow(ko_lfc) == 0L) {
  stop(
    "Knockout validation: no target genes found. ",
    "Check KO_TARGETS$gene_id against rownames(counts)."
  )
}

ko_lfc <- ko_lfc %>%
  mutate(
    mutant = recode(
      strain,
      "RB2060" = "argk-2",
      "RB2598" = "argk-4"
    ),
    target_title = recode(
      target_gene,
      "argk-2" = "italic('argk-2')~mRNA",
      "argk-4" = "italic('argk-4')~mRNA"
    ),
    stars = case_when(
      is.na(padj)  ~ NA_character_,
      padj < 1e-4  ~ "***",
      padj < 1e-3  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ NA_character_
    ),
    # Place significance marks just below negative bars, as in the thesis.
    star_y = ifelse(!is.na(stars) & log2FC < 0, log2FC - 0.16, log2FC + 0.16)
  )

cat("\nKnockout validation (all target x mutant comparisons):\n")
print(ko_lfc, n = Inf)

# Give the stars a little room below the most negative bar without changing data.
ko_ymin <- min(c(ko_lfc$log2FC, ko_lfc$star_y), na.rm = TRUE) - 0.12
ko_ymax <- max(0.12, max(ko_lfc$log2FC, na.rm = TRUE) + 0.12)

panelB <- ggplot(
  ko_lfc,
  aes(x = mutant, y = log2FC, fill = mutant)
) +
  geom_col(
    width = 0.56,
    colour = "black",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = stars, y = star_y),
    size = 4.1,
    fontface = "bold",
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ target_title,
    nrow = 1,
    labeller = label_parsed
  ) +
  scale_fill_manual(
    values = c(
      "argk-2" = GROUP_COLORS[["RB2060_untreated"]],
      "argk-4" = GROUP_COLORS[["RB2598_untreated"]]
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = c(
      "argk-2" = expression(italic("argk-2")^"-/-"),
      "argk-4" = expression(italic("argk-4")^"-/-")
    )
  ) +
  coord_cartesian(ylim = c(ko_ymin, ko_ymax), clip = "off") +
  labs(
    x = NULL,
    y = expression(Log[2]~FC)
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, colour = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.35),
    axis.text.x = element_text(size = 9.5, colour = "black"),
    axis.text.y = element_text(size = 9.5, colour = "black"),
    axis.title.y = element_text(size = 11, colour = "black"),
    panel.spacing.x = unit(0.55, "cm"),
    plot.margin = margin(5, 6, 5, 6)
  )

ggsave(
  file.path(OUT_DIR, "Figure1_PanelB_KO_validation.pdf"),
  panelB,
  width = 5.2,
  height = 3.3
)


# 5. Panel C -- Volcano plots 

contrast_grid <- expand.grid(
  strain    = c("RB2060", "RB2598"),
  treatment = c("untreated", "treated"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    pretty_strain = unname(GENOTYPE_LABELS[strain]),
    label    = paste0(pretty_strain, " vs N2 (", treatment, ")"),
    file_tag = paste0(strain, "_vs_N2_", treatment)
  )

make_volcano <- function(res, title_text, show_legend = FALSE) {
  df <- as.data.frame(res) %>%
    rownames_to_column("Gene") %>%
    mutate(
      sig = case_when(
        is.na(padj)                                       ~ "NA",
        padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Up",
        padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Down",
        TRUE                                              ~ "ns"
      ),
      sig = factor(sig, levels = c("Up", "Down", "ns", "NA"))
    )

  # Keep zero exactly centered on the x axis without modifying any values.
  finite_fc <- df$log2FoldChange[is.finite(df$log2FoldChange)]
  x_lim <- ceiling(max(abs(finite_fc), na.rm = TRUE))
  if (!is.finite(x_lim) || x_lim < 2) x_lim <- 2

  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) +
    geom_point(alpha = 0.65, size = 0.95, na.rm = TRUE) +
    scale_colour_manual(
      values = c("Up" = COL_UP, "Down" = COL_DOWN,
                 "ns" = "black", "NA" = COL_NA),
      labels = c("UP", "DN", "Not significant", "NA"),
      drop = FALSE
    ) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
               linetype = "dashed", linewidth = 0.35) +
    geom_hline(yintercept = -log10(PADJ_CUTOFF),
               linetype = "dashed", linewidth = 0.35) +
    coord_cartesian(xlim = c(-x_lim, x_lim)) +
    labs(
      title = title_text,
      x = expression(log[2]~FC),
      y = expression(-log[10]~"(adjusted p-value)"),
      colour = NULL
    ) +
    theme_classic(base_size = 10.5) +
    theme(
      plot.title = element_text(face = "bold", size = 10.5, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.4),
      legend.position = if (show_legend) c(0.82, 0.84) else "none",
      legend.background = element_rect(fill = alpha("white", 0.88), colour = "black", linewidth = 0.3),
      legend.key.height = unit(0.32, "cm"),
      legend.key.width = unit(0.32, "cm"),
      legend.text = element_text(size = 8.5),
      plot.margin = margin(4, 4, 4, 4)
    )
}

volcano_plots <- vector("list", nrow(contrast_grid))
for (i in seq_len(nrow(contrast_grid))) {
  row      <- contrast_grid[i, ]
  this_dds <- if (row$treatment == "untreated") dds_un else dds_tr
  res <- results(this_dds,
                 contrast = c("genotype", row$strain, "N2"),
                 alpha    = PADJ_CUTOFF)
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_symbol")
  write.csv(res_df,
            file = file.path(DEG_DIR,
                             paste0(row$strain, "_vs_N2_",
                                    row$treatment, "_DEGs.csv")),
            row.names = FALSE)
  
  p <- make_volcano(res, row$label, show_legend = (i == nrow(contrast_grid)))
  volcano_plots[[i]] <- p
  ggsave(file.path(OUT_DIR,
                   paste0("Figure1_PanelC_Volcano_", row$file_tag, ".pdf")),
         p, width = 5, height = 4.5)
}

res_n2 <- results(dds_n2,
                  contrast = c("treatment", "treated", "untreated"),
                  alpha    = PADJ_CUTOFF)
n2_df <- as.data.frame(res_n2) %>%
  rownames_to_column("gene_symbol")
write.csv(n2_df,
          file = file.path(DEG_DIR, "N2_treated_vs_untreated_DEGs.csv"),
          row.names = FALSE)

panelC <- wrap_plots(volcano_plots, ncol = 4)


#  6. Panel D - Treatment-response Venn diagram 

# DEG filter: padj < 0.05 AND |log2FC| > 1 (both thresholds match figure caption)
deg_filter <- function(df) {
  df %>%
    filter(!is.na(padj), !is.na(log2FoldChange),
           padj < PADJ_CUTOFF,
           abs(log2FoldChange) > LFC_CUTOFF) %>%
    pull(gene_symbol) %>%
    unique()
}
read_deg <- function(file) {
  read.csv(file.path(DEG_DIR, file)) %>%
    deg_filter()
}

#  Write within-mutant treatment DEG files (needed for treated Venn) 
# Same approach as dds_n2: fit a per-strain model, then extract treatment contrast.
make_within_strain_dds <- function(strain_label) {
  cd <- coldata[coldata$genotype == strain_label, , drop = FALSE]
  cd$treatment <- relevel(droplevels(cd$treatment), ref = "untreated")
  cc <- counts[, rownames(cd)]
  d  <- DESeqDataSetFromMatrix(countData = cc, colData = cd,
                               design = ~ treatment)
  d  <- d[rowSums(counts(d)) >= MIN_COUNT_SUM, ]
  DESeq(d)
}

for (strain in c("RB2060", "RB2598")) {
  dds_s <- make_within_strain_dds(strain)
  res_s <- results(dds_s,
                   contrast = c("treatment", "treated", "untreated"),
                   alpha    = PADJ_CUTOFF)
  df_s  <- as.data.frame(res_s) %>% rownames_to_column("gene_symbol")
  write.csv(df_s,
            file      = file.path(DEG_DIR,
                                  paste0(strain, "_treated_vs_untreated_DEGs.csv")),
            row.names = FALSE)
}

message("Re-extracting per-genotype treatment DEGs from full-factorial model with shrinkage...")
write_full_model_trt_deg <- function(out_file, coef = NULL, contrast_list = NULL) {
  res <- if (!is.null(coef)) {
    DESeq2::lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
  } else {
    DESeq2::lfcShrink(dds, contrast = contrast_list, type = "ashr", quiet = TRUE)
  }
  df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_symbol")
  write.csv(df, file = file.path(DEG_DIR, out_file), row.names = FALSE)
  n_sig <- sum(!is.na(df$padj) & df$padj < PADJ_CUTOFF &
                 abs(df$log2FoldChange) > LFC_CUTOFF)
  message(sprintf("  %-45s %5d DEGs (padj<%.2f, |LFC|>%.0f)",
                  out_file, n_sig, PADJ_CUTOFF, LFC_CUTOFF))
}
write_full_model_trt_deg("N2_treated_vs_untreated_DEGs.csv",
                         coef = "treatment_treated_vs_untreated")
write_full_model_trt_deg("RB2060_treated_vs_untreated_DEGs.csv",
                         contrast_list = list(c("treatment_treated_vs_untreated",
                                                "genotypeRB2060.treatmenttreated")))
write_full_model_trt_deg("RB2598_treated_vs_untreated_DEGs.csv",
                         contrast_list = list(c("treatment_treated_vs_untreated",
                                                "genotypeRB2598.treatmenttreated")))

# Load gene lists 
genes_argk2_un  <- read_deg("RB2060_vs_N2_untreated_DEGs.csv")      # argk-2 unt vs N2 unt
genes_argk4_un  <- read_deg("RB2598_vs_N2_untreated_DEGs.csv")      # argk-4 unt vs N2 unt
genes_n2_stress <- read_deg("N2_treated_vs_untreated_DEGs.csv")      # N2 tx vs N2 unt
genes_argk2_tx  <- read_deg("RB2060_treated_vs_untreated_DEGs.csv")  # argk-2 tx vs argk-2 unt
genes_argk4_tx  <- read_deg("RB2598_treated_vs_untreated_DEGs.csv")  # argk-4 tx vs argk-4 unt

cat("\nDEG set sizes (padj < 0.05, |log2FC| > 1):\n")
cat(sprintf("  argk-2(-/-) vs N2 untreated          : %d\n", length(genes_argk2_un)))
cat(sprintf("  argk-4(-/-) vs N2 untreated          : %d\n", length(genes_argk4_un)))
cat(sprintf("  N2 treated vs N2 untreated           : %d\n", length(genes_n2_stress)))
cat(sprintf("  argk-2(-/-) treated vs untreated     : %d\n", length(genes_argk2_tx)))
cat(sprintf("  argk-4(-/-) treated vs untreated     : %d\n", length(genes_argk4_tx)))

# Single 3-set Venn 
venn_tr <- list(
  "N2"          = genes_n2_stress,
  "argk-2(-/-)" = genes_argk2_tx,
  "argk-4(-/-)" = genes_argk4_tx
)

make_venn <- function(venn_list, low, high, title_text, subtitle_text = NULL) {
  # VennDiagram::venn.diagram returns a gList of grobs. wrap_elements(panel=)
  # turns it into a patchwork-compatible plot we can layout next to ggplots.
  n <- length(venn_list)
  fills <- if (n == 2) c(low, high) else c(low, high, "grey75")
  vd <- venn.diagram(
    x          = venn_list,
    filename   = NULL,
    fill       = fills,
    alpha      = 0.55,
    lwd        = 1,
    col        = "black",
    fontfamily      = "sans",
    cat.fontfamily  = "sans",
    cex             = 1.05,
    cat.cex         = 0.95,
    cat.default.pos = "outer",
    main            = title_text,
    main.cex        = 1.05,
    main.fontfamily = "sans",
    sub             = subtitle_text,
    sub.cex         = 0.75,
    sub.fontfamily  = "sans",
    margin          = 0.04,
    disable.logging = TRUE
  )
  # Clean up the per-call VennDiagram .log file that's silently created
  unlink(list.files(pattern = "^VennDiagram.*\\.log$"))
  patchwork::wrap_elements(panel = grid::gTree(children = vd))
}

p_venn_tr <- make_venn(
  venn_tr, COL_TR["low"], COL_TR["high"],
  title_text    = "Treatment response across genotypes",
  subtitle_text = "treated vs untreated within each genotype | padj < 0.05, |log\u2082FC| > 1"
)

ggsave(file.path(OUT_DIR, "Figure1_PanelD_Venn_Treatment.pdf"),
       p_venn_tr, width = 5.2, height = 3.4)

panelD <- p_venn_tr


# 7. Composite figure 
# Final journal layout:
#   A/B on the first row
#   C = four volcano plots in one row
#   D = compact treatment-response Venn centered beneath

panelA_tag <- panelA +
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 15))

panelB_tag <- panelB +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold", size = 15))

# panelC is already a four-plot patchwork. Add its label without converting it
# to another grob first.
panelC_tag <- panelC +
  patchwork::plot_annotation(
    tag_levels = NULL,
    theme = theme(plot.margin = margin(2, 2, 2, 2))
  )

# Add the C label as a narrow text column so the volcano patchwork remains intact.
c_label <- ggplot() +
  annotate("text", x = 0, y = 1, label = "C", fontface = "bold", size = 5.3,
           hjust = 0, vjust = 1) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

panelC_row <- c_label | panelC_tag
panelC_row <- panelC_row + plot_layout(widths = c(0.055, 0.945))

d_label <- ggplot() +
  annotate("text", x = 0, y = 1, label = "D", fontface = "bold", size = 5.3,
           hjust = 0, vjust = 1) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

venn_center <- plot_spacer() | panelD | plot_spacer()
venn_center <- venn_center + plot_layout(widths = c(0.55, 1.35, 0.55))

panelD_row <- d_label | venn_center
panelD_row <- panelD_row + plot_layout(widths = c(0.055, 0.945))

figure1 <- (panelA_tag | panelB_tag) /
  panelC_row /
  panelD_row +
  plot_layout(
    heights = c(1.0, 0.92, 0.60)
  )

ggsave(
  file.path(OUT_DIR, "Figure1_Composite.pdf"),
  figure1,
  width = 12,
  height = 8.9,
  bg = "white"
)

ggsave(
  file.path(OUT_DIR, "Figure1_Composite.png"),
  figure1,
  width = 12,
  height = 8.9,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(OUT_DIR, "Figure1_Composite_600dpi.tiff"),
  figure1,
  width = 12,
  height = 8.9,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

message("Done. All Figure 1 outputs written to: ", normalizePath(OUT_DIR))
