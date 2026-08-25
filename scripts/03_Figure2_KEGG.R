#!/usr/bin/env Rscript
# 03_Figure2_KEGG.R
#
# PURPOSE
#   Render Figure 2 -- KEGG functional composition of the CoCl2 stress response.
#   Panels A/B summarize three biologically interpretable main contrasts; Panel C
#   summarizes the genotype x treatment interaction term.
#
# INPUTS (all produced by 01_DESeq2_and_TF_enrichment.R, in outputs/DESeq2/)
#   DESeq2_treatment_in_N2.csv         -- N2: treated vs untreated
#   DESeq2_treatment_in_RB2060.csv     -- argk-2: treated vs untreated
#   DESeq2_genotype_untreated.csv      -- argk-2 vs N2 (untreated baseline)
#   DESeq2_GxE_interaction.csv         -- genotype x treatment interaction term
#
# Plus the offline KEGG cel mappings (in data/reference/):
#   kegg_cel_pathway_to_gene.csv
#   kegg_cel_pathway_names.csv
#
# DEG-calling filter for KEGG composition throughout: padj < 0.05 and |log2FC| > 1.
# The same threshold is applied consistently to Panels A/B and the G x E genes in Panel C.
#
# OUTPUTS
#   outputs/figures/Figure2/Figure2_KEGG_journal_ready.pdf   -- 3-panel figure
#   outputs/figures/Figure2/Figure2_KEGG_journal_ready.png   -- raster preview
#   outputs/DESeq2/KEGG_Functional_Composition_All_Contrasts.csv
#                                                        -- long-format DEG x pathway

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(stringr); library(purrr); library(ggplot2)
  library(patchwork)
})
for (utils_path in c("scripts/_utils.R", "_utils.R", "../scripts/_utils.R")) {
  if (file.exists(utils_path)) { source(utils_path); break }
}
stopifnot(exists("load_kegg_offline"))

DESEQ_DIR <- "outputs/DESeq2"
FIG_DIR   <- "outputs/figures/Figure2"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# DEG threshold used throughout Figure 2
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1
message(sprintf("Figure 2 DEG filter: padj < %.2f AND |log2FC| > %.1f",
                PADJ_CUTOFF, LFC_CUTOFF))

#  1. Load contrast tables and label them with friendly contrast names 
read_contrast <- function(filename, contrast_label) {
  path <- file.path(DESEQ_DIR, filename)
  if (!file.exists(path)) {
    stop("Missing ", path, "\n  Run 01_DESeq2_and_TF_enrichment.R first.")
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  # First column is the unnamed rownames col -> gene_symbol fallback
  if (colnames(df)[1] == "X" || colnames(df)[1] == "") {
    colnames(df)[1] <- "gene_symbol_rn"
  }
  if (!"gene_symbol" %in% colnames(df)) df$gene_symbol <- df$gene_symbol_rn
  df$contrast <- contrast_label
  df
}

deg_n2_trt       <- read_contrast("DESeq2_treatment_in_N2.csv",     "N2: Treated vs Untreated")
deg_rb2060_trt   <- read_contrast("DESeq2_treatment_in_RB2060.csv", "argk-2: Treated vs Untreated")
deg_geno_unt     <- read_contrast("DESeq2_genotype_untreated.csv",  "argk-2 vs N2 (Untreated)")
deg_interaction  <- read_contrast("DESeq2_GxE_interaction.csv",     "Interaction (G x E)")

message(sprintf("Loaded contrasts (rows): N2=%d, argk-2=%d, geno_unt=%d, int=%d",
                nrow(deg_n2_trt), nrow(deg_rb2060_trt),
                nrow(deg_geno_unt), nrow(deg_interaction)))

# 2. Build KEGG gene -> pathway map 
kegg_paths <- load_kegg_offline()
pathway_names <- as_tibble(kegg_paths$KEGGPATHID2NAME) |>
  setNames(c("pathway_id", "pathway_name"))
gene2path <- as_tibble(kegg_paths$KEGGPATHID2EXTID) |>
  setNames(c("pathway_id", "kegg_gene")) |>
  inner_join(pathway_names, by = "pathway_id") |>
  mutate(kegg_gene = str_remove(kegg_gene, "^CELE_"))
message(sprintf("KEGG: %d unique pathways, %d gene-pathway edges",
                nrow(pathway_names), nrow(gene2path)))

# 3. Compose DEG x KEGG cross-table 
# KEGG figure threshold: padj < 0.05 and |log2FC| > 1
make_deg_kegg <- function(df) {
  df |>
    transmute(Gene = gene_symbol,
              log2FC = log2FoldChange,
              padj   = padj,
              contrast = contrast) |>
    filter(!is.na(log2FC), !is.na(padj),
           padj < PADJ_CUTOFF, abs(log2FC) > LFC_CUTOFF) |>
    inner_join(gene2path, by = c("Gene" = "kegg_gene"),
               relationship = "many-to-many")
}

deg_kegg_3 <- bind_rows(
  make_deg_kegg(deg_n2_trt),
  make_deg_kegg(deg_rb2060_trt),
  make_deg_kegg(deg_geno_unt)
)
message(sprintf("DEG-to-KEGG rows across the three Panel A/B contrasts: %d",
                nrow(deg_kegg_3)))

# Save canonical KEGG composition CSV used by Panels A/B
write.csv(deg_kegg_3,
          file.path(DESEQ_DIR, "KEGG_Functional_Composition_All_Contrasts.csv"),
          row.names = FALSE)
message("Wrote: ", file.path(DESEQ_DIR, "KEGG_Functional_Composition_All_Contrasts.csv"))

# Reusable shortcut for Panel A code below
deg_kegg <- deg_kegg_3

# 
# 4. Panel A: KEGG violin -- THESIS-STYLE HORIZONTAL LAYOUT 
# 

message("Building Panel A in thesis-style layout...")

# Select pathways using the UPDATED data exactly as before.
pathway_order_A <- deg_kegg |>
  group_by(contrast, pathway_name) |>
  filter(n_distinct(Gene) >= 5) |>
  group_by(pathway_name) |>
  summarise(mean_abs = mean(abs(log2FC)), .groups = "drop") |>
  arrange(desc(mean_abs)) |>
  slice_head(n = 12) |>
  pull(pathway_name)


# Friendly labels chosen to reproduce the original thesis layout.
plot_A_data <- deg_kegg |>
  group_by(contrast, pathway_name) |>
  filter(n_distinct(Gene) >= 3) |>
  ungroup() |>
  filter(pathway_name %in% pathway_order_A) |>
  mutate(
    pathway_name = factor(
      pathway_name,
      levels = rev(pathway_order_A)
    ),
    
    contrast = case_when(
      contrast == "argk-2 vs N2 (Untreated)" ~
        "Genotype effect",
      
      contrast == "N2: Treated vs Untreated" ~
        "Treatment effect (N2)",
      
      contrast == "argk-2: Treated vs Untreated" ~
        "Treatment effect (argk-2−/−)",
      
      TRUE ~ contrast
    )
  )


# Preserve the same thesis color logic:
# blue = genotype effect
# red  = WT treatment effect
# green = argk-2 treatment effect
colors_A <- c(
  "Genotype effect"            = "#4E79A7",
  "Treatment effect (N2)"      = "#E15759",
  "Treatment effect (argk-2−/−)" = "#59A14F"
)


plot_A_data$contrast <- factor(
  plot_A_data$contrast,
  levels = names(colors_A)
)


# 
# IMPORTANT CHANGE:
# facet_grid(. ~ contrast) puts the three effects SIDE-BY-SIDE,
# reproducing the structure of the original thesis figure.
# 

panel_A <- ggplot(
  plot_A_data,
  aes(
    x = log2FC,
    y = pathway_name,
    fill = contrast
  )
) +
  
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.85,
    linewidth = 0.25
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.35
  ) +
  
  scale_fill_manual(
    values = colors_A
  ) +
  
  facet_grid(
    . ~ contrast,
    scales = "fixed",
    space = "fixed"
  ) +
  
  scale_x_continuous(
    breaks = c(-4, 0, 4)
  ) +
  
  labs(
    x = expression(log[2]~Fold~Change),
    y = NULL
  ) +
  
  theme_classic(base_size = 10) +
  
  theme(
    # pathway labels
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 9,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.title.y = element_blank(),
    
    # Thesis figure used simple titles rather than grey facet boxes.
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "plain",
      size = 10,
      color = "black"
    ),
    
    legend.position = "none",
    
    panel.spacing.x = unit(
      0.65,
      "cm"
    ),
    
    axis.line = element_line(
      linewidth = 0.35
    ),
    
    axis.ticks = element_line(
      linewidth = 0.3
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 3,
      l = 5
    )
  )



# 5. Panel B: DEG coverage heatmap (three main contrasts only)

message("Building Panel B...")

# Panel B deliberately excludes the G x E interaction. It is a compact summary
# of the same three padj < 0.05 and |log2FC| > 1 DEG sets shown in Panel A.
panelB_data <- deg_kegg_3 |>
  group_by(pathway_name, contrast) |>
  summarise(n_deg = n_distinct(Gene), .groups = "drop")

# Keep the 15 pathways with the greatest total DEG representation across the
# three contrasts. n_distinct(Gene) prevents duplicate KEGG edges from inflating
# the displayed number of DEGs.
top15_paths <- panelB_data |>
  group_by(pathway_name) |>
  summarise(total = sum(n_deg), .groups = "drop") |>
  arrange(desc(total), pathway_name) |>
  slice_head(n = 15) |>
  pull(pathway_name)

col_labels <- c(
  "Treat. eff. (N2)",
  "Treat. eff. (argk-2−/−)",
  "Genotype effect"
)

heatmap_df <- panelB_data |>
  filter(pathway_name %in% top15_paths) |>
  mutate(
    contrast_short = case_when(
      contrast == "N2: Treated vs Untreated" ~ "Treat. eff. (N2)",
      contrast == "argk-2: Treated vs Untreated" ~ "Treat. eff. (argk-2−/−)",
      contrast == "argk-2 vs N2 (Untreated)" ~ "Genotype effect",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(contrast_short)) |>
  complete(
    pathway_name = top15_paths,
    contrast_short = col_labels,
    fill = list(n_deg = 0)
  ) |>
  mutate(
    pathway_name = factor(pathway_name, levels = rev(top15_paths)),
    contrast_short = factor(contrast_short, levels = col_labels)
  )

write_csv(
  heatmap_df |> mutate(across(where(is.factor), as.character)),
  file.path(FIG_DIR, "Figure2_panelB_counts.csv")
)

panel_B <- ggplot(
  heatmap_df,
  aes(x = contrast_short, y = pathway_name, fill = n_deg)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = n_deg), size = 3.0, color = "black") +
  scale_fill_gradient(
    low = "#F0F4FA",
    high = "#1F4E79",
    name = "DEGs",
    limits = c(0, NA)
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 35, hjust = 1, vjust = 1, size = 9, color = "black"
    ),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.width = unit(0.35, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(3, 5, 3, 3)
  )


# 6. Panel C: Interaction-term DEGs 


message("Building Panel C...")

interaction_sig <- deg_interaction |>
  filter(
    !is.na(log2FoldChange),
    !is.na(padj),
    padj < PADJ_CUTOFF,
    abs(log2FoldChange) > LFC_CUTOFF
  ) |>
  transmute(
    Gene = gene_symbol,
    log2FC = log2FoldChange,
    padj
  )


interaction_kegg <- interaction_sig |>
  inner_join(
    gene2path,
    by = c("Gene" = "kegg_gene"),
    relationship = "many-to-many"
  )


interaction_kegg_filt <- interaction_kegg |>
  group_by(pathway_name) |>
  filter(n_distinct(Gene) >= 3) |>
  ungroup()


int_pathway_order <- interaction_kegg_filt |>
  group_by(pathway_name) |>
  summarise(
    n = n_distinct(Gene),
    .groups = "drop"
  ) |>
  arrange(desc(n)) |>
  slice_head(n = 15) |>
  pull(pathway_name)


bar_data_C <- interaction_kegg_filt |>
  filter(
    pathway_name %in% int_pathway_order
  ) |>
  
  mutate(
    direction = ifelse(
      log2FC > 0,
      "Positive G×E",
      "Negative G×E"
    ),
    
    pathway_name = factor(
      pathway_name,
      levels = rev(int_pathway_order)
    )
  ) |>
  
  group_by(
    pathway_name,
    direction
  ) |>
  
  summarise(
    n = n_distinct(Gene),
    .groups = "drop"
  )

write_csv(
  bar_data_C |> mutate(across(where(is.factor), as.character)),
  file.path(FIG_DIR, "Figure2_panelC_interaction_counts.csv")
)

panel_C <- ggplot(
  bar_data_C,
  aes(
    x = n,
    y = pathway_name,
    fill = direction
  )
) +
  
  geom_col(
    width = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      "Negative G×E" = "#2166AC",
      "Positive G×E" = "#B2182B"
    ),
    breaks = c(
      "Negative G×E",
      "Positive G×E"
    )
  ) +
  
  labs(
    x = "G×E interaction DEGs",
    y = NULL
  ) +
  
  theme_classic(base_size = 9) +
  
  theme(
    axis.text.y = element_text(
      size = 8.5,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 8.5,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 9.5,
      color = "black"
    ),
    
    axis.title.y = element_blank(),
    
    legend.position = "bottom",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.key.size = unit(
      0.3,
      "cm"
    ),
    
    legend.box.spacing = unit(
      0.1,
      "cm"
    ),
    
    plot.margin = margin(
      3, 5, 3, 3
    ),
    
    axis.line = element_line(
      linewidth = 0.35
    ),
    
    axis.ticks = element_line(
      linewidth = 0.3
    )
  )


# 7. Assemble Layout

bottom_row <- panel_B | panel_C


fig2 <- panel_A /
  bottom_row +
  
  plot_layout(
    heights = c(
      0.90,
      1.05
    ),
    
    widths = 1
  ) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    )
  )

out_pdf <- file.path(
  FIG_DIR,
  "Figure2_KEGG.pdf"
)

out_png <- file.path(
  FIG_DIR,
  "Figure2_KEGG.png"
)

out_tiff <- file.path(
  FIG_DIR,
  "Figure2_KEGG_600dpi.tiff"
)


ggsave(
  out_pdf,
  fig2,
  width = 10.5,
  height = 8.5,
  bg = "white"
)


ggsave(
  out_png,
  fig2,
  width = 10.5,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  out_tiff,
  fig2,
  width = 10.5,
  height = 8.5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

write_csv(
  interaction_kegg,
  file.path(FIG_DIR, "Figure2_panelC_interaction_gene_pathway_map.csv")
)


message(
  "\n✔ Figure 2 complete"
)
