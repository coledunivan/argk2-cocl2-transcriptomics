#!/usr/bin/env Rscript

# Fig3_combined.R
#
#
# PANEL A
#   GO BP enrichment: N2 treated vs untreated
#
# PANEL B
#   GO BP enrichment: argk-2-/- vs N2 in treated animals
#
# PANEL C
#   Immune effect map:
#   genotype effect vs treatment effect
#
# ANALYSIS DESIGN
#
# PANELS A/B:
#
#   Use the CANONICAL FULL-MODEL DESeq2 contrast tables produced by the main
#   DE pipeline:
#
#     outputs/DESeq2/DESeq2_treatment_in_N2.csv
#     outputs/DESeq2/DESeq2_genotype_treated.csv
#
#
#     CELE gene ID
#        ↓
#     celegans.db
#        ↓
#     ENTREZID
#        ↓
#     clusterProfiler::enrichGO()
#        ↓
#     org.Ce.eg.db
#
#   GO input:
#     DESeq2 padj < 0.10
#
#   GO:
#     BP
#     BH correction
#     pvalueCutoff = 0.10
#
#
# PANEL C:
#
#
#     ~ genotype * treatment +
#       genotype:experiment +
#       treatment:experiment
#
#   followed by the original dds_simple reparameterization:
#
#     ~ genotype + treatment + genotype:treatment
#
#   and the original:
#
#     TFLink → offline ENTREZ mapping
#     immune GO regex
#     effect-category classification
#
#   Only the single most-negative genotype-effect observation is omitted
#   FROM THE DISPLAY. It remains in Fig3_panel_C_data.csv.
#
#
# PLOTTING
#   * Thesis-style A | (B / C) geometry
#   * Only the long "biological process involved in interspecies interaction
#     between organisms" label is abbreviated
#   * Panel C symmetric around genotype effect = 0
#   * Panel C enlarged relative to previous recreation
#
#
# OUTPUTS
# outputs/figures/Figure3/
#
#   Fig3_FINAL_FULLMODEL_OLDGO_THESIS_STYLE.pdf
#   Fig3_FINAL_FULLMODEL_OLDGO_THESIS_STYLE.png
#   Fig3_panel_A_GO.csv
#   Fig3_panel_B_GO.csv
#   Fig3_panel_C_data.csv
#

# 0. Packages

suppressPackageStartupMessages({

  library(DESeq2)

  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)

  library(ggplot2)
  library(patchwork)
  library(scales)

  # Original thesis GO tools
  library(clusterProfiler)
  library(org.Ce.eg.db)
  library(AnnotationDbi)
  library(celegans.db)
  library(GO.db)
})

#  figure aesthetics 
# All typography, line weights, colour and output geometry resolve to
# scripts/_theme.R. Edit that file to restyle; nothing here changes data.
if (file.exists("scripts/_theme.R")) source("scripts/_theme.R") else source("_theme.R")

for (utils_path in c("scripts/_utils.R", "_utils.R", "../scripts/_utils.R")) {

  if (file.exists(utils_path)) {

    source(utils_path)

    break
  }
}

stopifnot(exists("load_go_offline"))

# Avoid namespace conflicts
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
group_by  <- dplyr::group_by
distinct  <- dplyr::distinct
summarise <- dplyr::summarise
summarize <- dplyr::summarize

# 1. Configuration

DESEQ_DIR <- "outputs/DESeq2"

OUT_DIR <- "outputs/figures/Figure3"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# GO threshold
GO_DE_PADJ <- 0.10
GO_P_CUT   <- 0.10

# Panel C significance threshold
DE_PADJ_CUT <- 0.05

MIN_SUM <- 10

PANEL_A_TOP <- 20

# 2. Locate Panel C source files

find_first <- function(patterns, label) {

  hits <- character(0)

  for (p in patterns) {

    hits <- c(hits, Sys.glob(p))
  }

  hits <- unique(hits)

  if (!length(hits)) {

    stop(sprintf("Could not find %s file. Tried: %s", label, paste(patterns, collapse = ", ")))
  }

  hits[1]
}

rna_file_C <- find_first(c("data/reference/RNAseq_to_TF_Targets*.csv",
  "data/reference/RNAseq_to_TF_Targets*.csv.gz", "RNAseq_to_TF_Targets*.csv",
  "RNAseq_to_TF_Targets*.csv.gz", "*TF_Targets*.csv"), "Panel C count/TFLink")

meta_file <- find_first(c("data/Sample_Metadata_Table*.csv", "Sample_Metadata_Table*.csv",
  "*Metadata*.csv"), "metadata")

message("Panel C source: ", rna_file_C)

message("Metadata:       ", meta_file)

# 3. Immune GO regex — ORIGINAL Panel C logic

immune_regex <- paste0("(immune|innate|host defen[cs]e|antimicrob|pathogen|",
  "bacteri|fungi?|virus|viral|", "antiviral|antibacteri|antifung|",
  "lysozyme|lys-\\d+|clec|c[- ]?type lectin|lectin|",
  "irg-\\d+|pals-\\d+|nlp-\\d+|spp-\\d+|defensin|", "pattern[- ]?recognition|pgrp|peptidoglycan|",
  "tlr|toll[- ]like|nlr|nod[- ]like|", "interferon|ifn-|complement|",
  "response to bacteri|response to fung|response to virus)")

# PART I
# PANELS A/B
# 4. Locate original RNA-seq count matrix

rna_file_GO <- find_first(c("data/RNASEQ61125.csv", "RNASEQ61125.csv",
  "data/reference/RNASEQ61125.csv"), "original GO count matrix")

message("GO counts: ", rna_file_GO)

# 5. Load original GO count matrix

raw_counts_GO <- read.csv(rna_file_GO, check.names = FALSE, stringsAsFactors = FALSE)

colnames(raw_counts_GO)[1] <- "Gene"

raw_counts_GO <- raw_counts_GO %>% filter(!duplicated(Gene))

counts_GO <- raw_counts_GO %>% column_to_rownames("Gene")

meta_GO <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)

meta_GO$SampleLabel <- trimws(meta_GO$SampleLabel)

rownames(meta_GO) <- meta_GO$SampleLabel

# 6. Align sample names

message("Aligning GO counts to metadata...")

shared_GO <- intersect(colnames(counts_GO), rownames(meta_GO))

if (length(shared_GO) == 0L) {

  message("No direct SampleLabel overlap. Trying alternate sample names...")

  if (!exists("harmonize_sample_names")) {

    for (utils_path in c("scripts/_utils.R", "_utils.R", "../scripts/_utils.R")) {

      if (file.exists(utils_path)) {

        source(utils_path)

        break
      }
    }
  }

  if (!exists("harmonize_sample_names")) {

    stop("Could not find harmonize_sample_names() in _utils.R.")
  }

  alt_col_GO <- intersect(c("GEO_library_name", "library_name", "Sample", "sample"),
    colnames(meta_GO))[1]

  if (is.na(alt_col_GO)) {

    stop("No alternate metadata sample-name column found.")
  }

  message("Using alternate sample-name column: ", alt_col_GO)

  counts_GO <- harmonize_sample_names(as.matrix(counts_GO), meta_GO, target_col = "SampleLabel",
    alt_col = alt_col_GO)

  shared_GO <- intersect(colnames(counts_GO), rownames(meta_GO))
}

if (length(shared_GO) == 0L) {

  stop("Still no overlapping GO samples after harmonization.")
}

message("GO shared samples: ", length(shared_GO))

counts_GO <- counts_GO[, shared_GO, drop = FALSE]

meta_GO <- meta_GO[shared_GO,, drop = FALSE]

counts_GO <- as.matrix(counts_GO)

storage.mode(counts_GO) <- "integer"

stopifnot(identical(colnames(counts_GO), rownames(meta_GO)))

#reference levels
# Treated is deliberately the treatment reference.
#
# Therefore:
#
# genotype_RB2060_vs_N2
#
# = RB2060 vs N2 specifically in treated animals.
#

meta_GO$genotype <- relevel(factor(meta_GO$genotype), ref = "N2")

meta_GO$treatment <- relevel(factor(meta_GO$treatment), ref = "treated")

meta_GO$experiment <- factor(meta_GO$experiment)

# 8. REDUCED thesis DESeq2 model
# Experiment is modeled as an additive nuisance/batch effect.

message("Running reduced GO model: ~ genotype * treatment + experiment")

dds_GO <- DESeqDataSetFromMatrix(countData = counts_GO, colData   = meta_GO,
  design    = ~ genotype * treatment + experiment)

dds_GO <- dds_GO[rowSums(counts(dds_GO)) > MIN_SUM,]

dds_GO <- DESeq(dds_GO)

message("GO model coefficient names:")

print(resultsNames(dds_GO))

# 9. Extract A/B coefficients exactly as thesis

if (!requireNamespace("apeglm", quietly = TRUE)) {

  stop("Package 'apeglm' is required.")
}

# Panel A
res_GO_A <- lfcShrink(dds_GO, coef = "treatment_untreated_vs_treated", type = "apeglm")

# Panel B

res_GO_B <- lfcShrink(dds_GO, coef = "genotype_RB2060_vs_N2", type = "apeglm")

# 10. CELE -> ENTREZ mapping

gene_df_GO <- data.frame(Gene = rownames(dds_GO), stringsAsFactors = FALSE)

gene_df_GO$Stripped <- sub("^CELE_", "", gene_df_GO$Gene)

entrez_map_GO <- suppressMessages(AnnotationDbi::select(celegans.db, keys = gene_df_GO$Stripped,
  columns = c("WORMBASE", "ENTREZID"), keytype = "SYMBOL"))

entrez_map_GO <- entrez_map_GO %>% mutate(Gene = paste0("CELE_", SYMBOL)) %>% select(Gene,
  WORMBASE, ENTREZID)

gene_df_GO <- left_join(gene_df_GO, entrez_map_GO, by = "Gene")

message(sprintf("GO ENTREZ mapping: %d / %d genes", sum(!is.na(gene_df_GO$ENTREZID)),
  nrow(gene_df_GO)))

#clusterProfiler GO enrichment

run_original_GO_BP <- function(res_obj, label) {

  message("\nGO enrichment for: ", label)

  res_df <- as.data.frame(res_obj) %>% rownames_to_column("Gene") %>% left_join(gene_df_GO,
    by = "Gene") %>% filter(!is.na(ENTREZID), !is.na(padj), padj < GO_DE_PADJ)

  entrez_ids <- unique(res_df$ENTREZID)

  message(sprintf("  Significant mapped genes: %d", length(entrez_ids)))

  if (length(entrez_ids) < 5) {

    warning("Not enough significant genes for ", label)

    return(NULL)
  }

  ego <- suppressMessages(clusterProfiler::enrichGO(gene          = entrez_ids,
    OrgDb         = org.Ce.eg.db, keyType       = "ENTREZID", ont           = "BP",
    pAdjustMethod = "BH", pvalueCutoff  = GO_P_CUT, readable      = TRUE))

  if (is.null(ego) || nrow(ego@result) == 0) {

    warning("No significant GO BP terms for ", label)

    return(NULL)
  }

  as.data.frame(ego) %>% mutate(GeneRatio_num = vapply(strsplit(GeneRatio, "/"), function(z) {

          as.numeric(z[1]) / as.numeric(z[2])
        }, numeric(1)))
}

# 12. Run A/B enrichment

go_A <- run_original_GO_BP(res_GO_A, "Panel A: N2 treated vs untreated")

go_B <- run_original_GO_BP(res_GO_B, "Panel B: argk-2-/- vs N2 (treated)")

if (!is.null(go_A)) {

  write.csv(go_A, file.path(OUT_DIR, "Fig3_panel_A_GO.csv"), row.names = FALSE)
}

if (!is.null(go_B)) {

  write.csv(go_B, file.path(OUT_DIR, "Fig3_panel_B_GO.csv"), row.names = FALSE)
}

# 13. Panel B thesis immune hierarchy

PANEL_B_TERMS <- c("immune response", "immune system process",
  "defense response to other organism", "response to biotic stimulus",
  "response to external biotic stimulus", "response to other organism", "defense response",
  "innate immune response", "defense response to symbiont")

# PART II
# PANEL C

# 10. Load original Panel C source
rna_raw_C <- read.csv(rna_file_C, check.names = FALSE, stringsAsFactors = FALSE)

meta_C <- read.csv(meta_file, stringsAsFactors = FALSE)

count_cols_C <- grep("^A\\d+", colnames(rna_raw_C), value = TRUE)

stopifnot(length(count_cols_C) > 0)

counts_C <- rna_raw_C %>% select(Gene, all_of(count_cols_C)) %>% distinct() %>%
  column_to_rownames("Gene")

meta_C <- meta_C %>% filter(SampleLabel %in% colnames(counts_C)) %>% column_to_rownames("SampleLabel"
  ) %>% mutate(genotype = factor(genotype, levels = c("N2", "RB2060", "RB2598")),
  treatment = factor(treatment, levels = c("untreated", "treated"))) %>% filter(!is.na(genotype),
  !is.na(treatment), !is.na(experiment))

shared_C <- intersect(colnames(counts_C), rownames(meta_C))

if (!length(shared_C)) {

  stop("No overlapping samples for Panel C.")
}

counts_C <- counts_C[, shared_C, drop = FALSE]

meta_C <- meta_C[shared_C,, drop = FALSE]

# 11. Panel C DESeq2 workflow

message("Running original Panel C DESeq2 computation...")

dds_C <- DESeqDataSetFromMatrix(countData = counts_C, colData = meta_C, design = ~ genotype * treatment +
  genotype:experiment + treatment:experiment)

dds_C <- dds_C[rowSums(counts(dds_C)) > MIN_SUM,]

dds_C <- DESeq(dds_C)

# Preserve the original second DESeq2 model exactly.

dds_simple_C <- dds_C

dds_simple_C$genotype <- relevel(dds_simple_C$genotype, ref = "N2")

dds_simple_C$treatment <- relevel(dds_simple_C$treatment, ref = "untreated")

design(dds_simple_C) <- ~ genotype + treatment + genotype:treatment

dds_simple_C <- DESeq(dds_simple_C)

# 12. Panel C contrast helpers

get_genotype_within_treatment <- function(dds_obj, genotype_alt = "RB2060", trt = "untreated") {

  rn <- resultsNames(dds_obj)

  main_coef <- rn[grepl("^genotype", rn) & grepl(genotype_alt, rn) & !grepl("treatment", rn)]

  inter <- rn[grepl("genotype", rn) & grepl(genotype_alt, rn) & grepl("treatment.*treated", rn)]

  if (trt == "untreated") {

    results(dds_obj, name = main_coef)

  } else {

    results(dds_obj, list(c(main_coef, inter)))
  }
}

get_treatment_within_genotype <- function(dds_obj, genotype_level = "N2") {

  rn <- resultsNames(dds_obj)

  trt_pos <- rn[grepl("^treatment", rn) & grepl("treated", rn)]

  inter <- rn[grepl("genotype", rn) & grepl(genotype_level, rn) & grepl("treatment.*treated", rn)]

  if (genotype_level == "N2") {

    results(dds_obj, name = trt_pos)

  } else {

    results(dds_obj, list(c(trt_pos, inter)))
  }
}

# 13. Extract Panel C contrasts

message("Extracting original Panel C contrasts...")

res_treat_N2_C <- get_treatment_within_genotype(dds_simple_C, "N2")

res_treat_RB_C <- get_treatment_within_genotype(dds_simple_C, "RB2060")

res_geno_unt_C <- get_genotype_within_treatment(dds_simple_C, "RB2060", "untreated")

res_geno_trt_C <- get_genotype_within_treatment(dds_simple_C, "RB2060", "treated")

inter_name_C <- grep("genotype.*RB2060.*treatment.*treated", resultsNames(dds_simple_C),
  value = TRUE)

res_gxe_C <- results(dds_simple_C, name = inter_name_C)

# 14 Panel C offline gene → ENTREZ mapping

message("Running original Panel C TFLink/offline ENTREZ mapping...")

go_map_full_C <- load_go_offline("data/reference/org_Ce_eg_GO_map.csv")

sym_to_entrez_map_C <- go_map_full_C %>% distinct(SYMBOL, ENTREZID) %>% filter(!is.na(SYMBOL),
  !is.na(ENTREZID))

sym_to_entrez_map_C$SYMBOL_lc <- tolower(sym_to_entrez_map_C$SYMBOL)

tfl_raw_C <- read.csv(resolve_tflink_path(), stringsAsFactors = FALSE, check.names = FALSE)

tfl_raw_C <- tfl_raw_C[colnames(tfl_raw_C) != "" & !is.na(colnames(tfl_raw_C)),, drop = FALSE]

cosmid_to_sym_C <- data.frame(gene = trimws(as.character(tfl_raw_C$Gene)), SYMBOL_lc = tolower(trimws(as.character(tfl_raw_C$Name.Target
  ))), stringsAsFactors = FALSE)

cosmid_to_sym_C <- cosmid_to_sym_C[cosmid_to_sym_C$SYMBOL_lc != "" &
  cosmid_to_sym_C$SYMBOL_lc != "-" & !is.na(cosmid_to_sym_C$SYMBOL_lc),, drop = FALSE]

cosmid_to_sym_C <- cosmid_to_sym_C[!duplicated(cosmid_to_sym_C$gene),]

map_to_entrez_C <- function(gene_names) {

  cleaned <- sub("^CELE_", "", as.character(gene_names))

  step1 <- data.frame(gene = cleaned, stringsAsFactors = FALSE)

  step1 <- merge(step1, cosmid_to_sym_C, by = "gene", all.x = TRUE, sort = FALSE)

  step2 <- merge(step1, sym_to_entrez_map_C[, c("SYMBOL_lc", "ENTREZID")], by = "SYMBOL_lc",
    all.x = TRUE, sort = FALSE)

  needs_fallback <- is.na(step2$ENTREZID)

  if (any(needs_fallback)) {

    fb_keys <- tolower(step2$gene[needs_fallback])

    fb <- sym_to_entrez_map_C$ENTREZID[match(fb_keys, sym_to_entrez_map_C$SYMBOL_lc)]

    step2$ENTREZID[needs_fallback] <- fb
  }

  result <- step2$ENTREZID[match(cleaned, step2$gene)]

  names(result) <- gene_names

  result
}

all_genes_C <- rownames(counts_C)

gene_to_entrez_C <- map_to_entrez_C(all_genes_C)

message(sprintf("Panel C ENTREZ mapping: %d / %d genes", sum(!is.na(gene_to_entrez_C)), length(all_genes_C
  )))

# 15. Panel C master table

master_C <- data.frame(Gene = rownames(res_treat_N2_C), lfc_treat_N2 = res_treat_N2_C$log2FoldChange,
  lfc_treat_RB = res_treat_RB_C$log2FoldChange, lfc_geno_unt = res_geno_unt_C$log2FoldChange,
  lfc_geno_trt = res_geno_trt_C$log2FoldChange, lfc_gxe = res_gxe_C$log2FoldChange,
  padj_treat_N2 = res_treat_N2_C$padj, padj_treat_RB = res_treat_RB_C$padj, padj_geno_unt = res_geno_unt_C$padj,
  padj_geno_trt = res_geno_trt_C$padj, padj_gxe = res_gxe_C$padj, stringsAsFactors = FALSE)

# 16. Panel C classification

classified_C <- master_C %>% mutate(sig_treat = (!is.na(padj_treat_N2) & padj_treat_N2 <
  DE_PADJ_CUT) | (!is.na(padj_treat_RB) & padj_treat_RB < DE_PADJ_CUT), sig_geno = (!is.na(padj_geno_unt
  ) & padj_geno_unt < DE_PADJ_CUT) | (!is.na(padj_geno_trt) & padj_geno_trt < DE_PADJ_CUT),
  sig_inter = !is.na(padj_gxe) & padj_gxe < DE_PADJ_CUT, category = case_when(sig_inter | (sig_geno &
  sig_treat) ~ "Interaction", sig_geno & !sig_treat ~ "Genotype only", sig_treat & !sig_geno ~ "Treatment only",
  TRUE ~ NA_character_), mean_geno_lfc = rowMeans(cbind(lfc_geno_unt, lfc_geno_trt), na.rm = TRUE
  ), mean_treat_lfc = rowMeans(cbind(lfc_treat_N2, lfc_treat_RB), na.rm = TRUE)) %>% filter(!is.na(category
  ))

# 17. Panel C immune GO annotation

classified_C$ENTREZ <- gene_to_entrez_C[classified_C$Gene]

entrez_classified_C <- na.omit(unique(classified_C$ENTREZ))

go_tbl_C <- go_map_full_C %>% filter(ENTREZID %in% entrez_classified_C, ONTOLOGY == "BP",
  !is.na(GO))

term_tbl_C <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(go_tbl_C$GO),
  columns = "TERM", keytype = "GOID"))

go_annot_C <- go_tbl_C %>% left_join(term_tbl_C, by = c("GO" = "GOID")) %>% filter(!is.na(TERM))

immune_entrez_C <- go_annot_C %>% filter(str_detect(TERM, regex(immune_regex, ignore_case = TRUE))
  ) %>% pull(ENTREZID) %>% unique()

immune_df <- classified_C %>% filter(ENTREZ %in% immune_entrez_C) %>% distinct(Gene,
  .keep_all = TRUE)

message(sprintf("Immune-annotated genes for Panel C: %d", nrow(immune_df)))

message(sprintf("  by category: %s", paste(names(table(immune_df$category)), "=", table(immune_df$category
  ), collapse = ", ")))

# Save BEFORE display-only outlier removal
write.csv(immune_df, file.path(OUT_DIR, "Fig3_panel_C_data.csv"), row.names = FALSE)

# PART III
# PLOTTING
# 18.  GO dotplot

# Only one label is abbreviated:
# "biological process involved in interspecies interaction between organisms"
#           ↓
# "interspecies interaction"

make_go_dotplot <- function(go_df, title, top_n, panel_type = c("A", "B")) {

  panel_type <- match.arg(panel_type)

  if (is.null(go_df) || nrow(go_df) == 0) {

    return(ggplot() + theme_void() + labs(title = "No significant GO terms"))
  }

  # ---------------------------------------------------------------------------
  # Select displayed terms
  # ---------------------------------------------------------------------------

  if (panel_type == "B") {

    d <- go_df %>% filter(Description %in% PANEL_B_TERMS) %>% mutate(thesis_order = match(Description,
      PANEL_B_TERMS)) %>% arrange(thesis_order)

  } else {

    # Panel A keeps top 20 by GeneRatio
    d <- go_df %>% arrange(desc(GeneRatio_num)) %>% slice_head(n = top_n)
  }

  # Display-only GO abbreviation

  d <- d %>% mutate(Description_plot = case_when(Description == "biological process involved in interspecies interaction between organisms" ~ "interspecies interaction",
    Description == "biological process involved in interspecies interaction" ~ "interspecies interaction",
    TRUE ~ Description))

  d$Description_plot <- factor(d$Description_plot, levels = rev(d$Description_plot))

  # Plot

  p <- ggplot(d, aes(x = GeneRatio_num, y = Description_plot, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.95) + scale_color_gradient(low = "#E76F61", high = "#4C78A8",
    name = "p.adjust", labels = function(x) {

        format(x, scientific = TRUE, digits = 2)
      }, guide = guide_colorbar(barheight = unit(1.3, "cm"), barwidth = unit(0.28, "cm"))) + labs(x = "GeneRatio",
        y = NULL, title = title) + theme_bw(base_size = 10) + theme(panel.border = element_rect(colour = "black",
        fill = NA, linewidth = G3$lw_geom), panel.grid.major = element_line(colour = "grey90",
        linewidth = G3$lw_axis), panel.grid.minor = element_blank(), plot.title = element_text(size = 11,
        hjust = 0.5, face = "plain", margin = margin(b = 4)), axis.text.y = element_text(size = ifelse(panel_type == "A",
        8.4, 8.5), colour = "black", lineheight = 0.92), axis.text.x = element_text(size = g3_pt(9),
        colour = "black"), axis.title.x = element_text(size = 10, colour = "black",
        margin = margin(t = 3)), axis.ticks = element_line(colour = "black",
        linewidth = G3$lw_geom), legend.position = "right", legend.title = element_text(size = g3_pt(8.5)
        ), legend.text = element_text(size = g3_pt(8)), legend.key.size = unit(0.30, "cm"),
        legend.spacing.y = unit(0.02, "cm"), plot.margin = margin(3, 3, 3, 3))

  if (panel_type == "A") {

    p <- p + scale_size_continuous(range = c(2.7, 6.7), name = "Count")

  } else {

    p <- p + scale_size_continuous(range = c(2.6, 6.1), name = "Count",
      breaks = function(l) unique(floor(pretty(l)))) +
      scale_x_continuous(n.breaks = 3)
  }

  p
}

panel_A <- make_go_dotplot(go_A, "N2 Treated vs. Untreated", PANEL_A_TOP, panel_type = "A")

panel_B <- make_go_dotplot(go_B, expression(italic("argk-2")^"-/-" * " vs. N2 (Treated)"), length(PANEL_B_TERMS
  ), panel_type = "B")

pal_effect <- c("Genotype only" = "#009E73", "Treatment only" = "#0072B2", "Interaction" = "#D55E00"
  )

if (nrow(immune_df) < 2) {

  stop("Panel C original computation produced fewer than 2 immune genes.")
}

# Exactly ONE observation:
# the most-negative genotype effect.
outlier_idx <- which.min(immune_df$mean_geno_lfc)

outlier_gene <- immune_df$Gene[outlier_idx]

outlier_value <- immune_df$mean_geno_lfc[outlier_idx]

message(sprintf(paste0("Panel C display-only outlier removed: ", "%s (genotype effect = %.3f)"),
  outlier_gene, outlier_value))

immune_df_plot <- immune_df[-outlier_idx,, drop = FALSE]

remaining_max_abs <- max(abs(immune_df_plot$mean_geno_lfc), na.rm = TRUE)

# Clean half-unit ceiling
x_lim <- ceiling(remaining_max_abs * 2) / 2

# Maintain roughly thesis-like minimum width
x_lim <- max(1.5, x_lim)

message(sprintf("Panel C centered x-axis: %.1f to %.1f", -x_lim, x_lim))

x_breaks <- pretty(c(-x_lim, x_lim), n = 5)

x_breaks <- x_breaks[x_breaks >= -x_lim & x_breaks <= x_lim]

panel_C <- ggplot(immune_df_plot, aes(x = mean_geno_lfc, y = mean_treat_lfc, color = category)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = G3$lw_axis) + geom_vline(xintercept = 0,
  colour = "grey55", linewidth = G3$lw_axis) + geom_point(size = 2.3, alpha = 0.75, stroke = 0) +
  scale_color_manual(values = pal_effect, breaks = c("Genotype only", "Treatment only",
  "Interaction"), name = NULL) + scale_x_continuous(limits = c(-x_lim, x_lim), breaks = x_breaks,
  expand = expansion(mult = c(0.02, 0.02))) +
  # Thesis-like y display.
  # coord_cartesian avoids deleting underlying rows.
  coord_cartesian(ylim = c(-3, 5)) + scale_y_continuous(breaks = c(-2.5, 0, 2.5, 5)) + labs(x = expression("Genotype effect (Log"[2] *
    "FC)"), y = expression("Treatment effect (Log"[2] * "FC)")) + theme_bw(base_size = 10) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = G3$lw_geom),
    panel.grid.major = element_line(colour = "grey90", linewidth = G3$lw_axis), panel.grid.minor = element_blank(),
    axis.text = element_text(size = g3_pt(9), colour = "black"), axis.title.x = element_text(size = 10,
    colour = "black", margin = margin(t = 4)), axis.title.y = element_text(size = 10,
    colour = "black", margin = margin(r = 4)), axis.ticks = element_line(linewidth = G3$lw_geom,
    colour = "black"), legend.position = "right", legend.text = element_text(size = g3_pt(9)),
    legend.key.size = unit(0.28, "cm"), legend.spacing.y = unit(0.02, "cm"),
    # Smaller margins so graph occupies more of its allocated panel.
    plot.margin = margin(2, 2, 2, 2))

design <- "
AB
AC
"

final <- panel_A + panel_B + panel_C + plot_layout(design = design, widths = c(1.45, 1.10),
  heights = c(0.46, 0.54)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold",
  size = g3_pt(18), colour = "black"))

# 24. Save

out_pdf <- file.path(OUT_DIR, "Fig3.pdf")

out_png <- file.path(OUT_DIR, "Fig3.png")

ggsave(out_pdf, final, width = G3$w_double, height = 5.6, device = cairo_pdf, bg = "white")

ggsave(out_png, final, width = G3$w_double, height = 5.6, dpi = 500, bg = "white")

# 25. Final diagnostics

message("\n============================================================")

message("FIGURE 3 COMPLETE")

message("============================================================")

message("\nPanels A/B:")

message("  DE source = canonical FULL-MODEL DESeq2 outputs")

message("  GO = original thesis clusterProfiler::enrichGO")

message("  mapping = celegans.db -> ENTREZID")

message("  GO ontology = BP")

message("  GO adjustment = BH")

message("  GO DEG input = DESeq2 padj < ", GO_DE_PADJ)

message("\nPanel A GO terms returned: ", ifelse(is.null(go_A), 0, nrow(go_A)))

message("Panel B GO terms returned before display filtering: ", ifelse(is.null(go_B), 0,
  nrow(go_B)))

if (!is.null(go_B)) {

  n_B_display <- sum(go_B$Description %in% PANEL_B_TERMS)

  message("Panel B thesis immune terms displayed: ", n_B_display, " / ", length(PANEL_B_TERMS))

  missing_B_terms <- setdiff(PANEL_B_TERMS, go_B$Description)

  if (length(missing_B_terms) > 0) {

    message("Panel B thesis terms not significant/present under full model:")

    message("  ", paste(missing_B_terms, collapse = "\n  "))
  }
}

message("\nPanel C:")

message("  original recreated Figure 3 computation preserved")

message("  immune genes = ", nrow(immune_df))

message("  plotted genes = ", nrow(immune_df_plot))

message("  display-only excluded gene = ", outlier_gene)

message(sprintf("  centered x-axis = %.1f to %.1f", -x_lim, x_lim))

message("\nSaved PDF: ", out_pdf)

message("Saved PNG: ", out_png)
