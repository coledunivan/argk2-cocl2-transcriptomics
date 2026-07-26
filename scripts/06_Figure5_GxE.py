#!/usr/bin/env python3
"""
Figure 5 — minimal reproducible ARGK-2 G×E pipeline.

Canonical repo inputs:
    data/reference/RNAseq_to_TF_Targets.csv
    data/Sample_Metadata_Table.csv

This script:
  1. rebuilds the Figure 5 DESeq2 effect table;
  2. recreates the original Figure 5 panels A/B/C.

Deliberately omitted:
  STRING/PPI queries, nodes/edges, network metrics, TF enrichment tests,
  GRN summaries, alternate metadata files, and unrelated outputs.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from math import comb
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_PATH = Path(__file__).resolve()
ROOT = SCRIPT_PATH.parents[1] if SCRIPT_PATH.parent.name == "scripts" else Path.cwd()

def _resolve_tflink(root):
    """Prefer the plain CSV; fall back to the gzipped copy shipped on GitHub."""
    for name in ("RNAseq_to_TF_Targets.csv", "RNAseq_to_TF_Targets.csv.gz"):
        p = root / "data/reference" / name
        if p.exists():
            return p
    raise FileNotFoundError("data/reference/RNAseq_to_TF_Targets.csv[.gz] not found")


RNA_FILE  = _resolve_tflink(ROOT)
META_FILE = ROOT / "data/Sample_Metadata_Table.csv"

DESEQ_DIR = ROOT / "outputs/DESeq2"
FIG_DIR   = ROOT / "outputs/figures/Figure5"

DESEQ_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

EFFECT_FILE = DESEQ_DIR / "DESeq2_effect_classification.csv"
OUT_PDF = FIG_DIR / "Figure5_ARGK2_GxE.pdf"
OUT_PNG = FIG_DIR / "Figure5_ARGK2_GxE.png"

for required in (RNA_FILE, META_FILE):
    if not required.exists():
        sys.exit(f"[error] Missing required repo input: {required}")

R_EFFECT_SCRIPT = '\nsuppressPackageStartupMessages({\n  library(DESeq2)\n  library(dplyr)\n})\n\nRNA_FILE  <- "RNAseq_to_TF_Targets.csv"\nMETA_FILE <- "Sample_Metadata_Table.csv"\n\ncat("Loading canonical Figure 5 inputs...\\n")\n\nmerged <- read.csv(\n  RNA_FILE,\n  stringsAsFactors = FALSE,\n  check.names = FALSE\n)\n\ncolnames(merged) <- trimws(colnames(merged))\nmerged <- merged[, colnames(merged) != "", drop = FALSE]\n\nsample_cols <- grep("^A[0-9]+$", colnames(merged), value = TRUE)\n\nif (length(sample_cols) == 0L) {\n  stop("No A1..An sample columns found in RNAseq_to_TF_Targets.csv")\n}\n\ncounts_dedup <- merged[!duplicated(merged$Gene), , drop = FALSE]\ncount_mat_raw <- as.matrix(counts_dedup[, sample_cols, drop = FALSE])\nrownames(count_mat_raw) <- counts_dedup$Gene\nstorage.mode(count_mat_raw) <- "double"\n\nmeta <- read.csv(\n  META_FILE,\n  stringsAsFactors = FALSE,\n  check.names = FALSE\n)\n\ncolnames(meta) <- trimws(colnames(meta))\nmeta$SampleLabel <- trimws(meta$SampleLabel)\n\nmeta_sub <- meta[\n  meta$genotype %in% c("N2", "RB2060"),\n  ,\n  drop = FALSE\n]\n\nmeta_sub <- droplevels(meta_sub)\nrownames(meta_sub) <- meta_sub$SampleLabel\n\nmissing_in_counts <- setdiff(meta_sub$SampleLabel, colnames(count_mat_raw))\n\nif (length(missing_in_counts) > 0L) {\n  stop(\n    "Metadata samples missing from count matrix: ",\n    paste(missing_in_counts, collapse = ", ")\n  )\n}\n\ncount_mat <- count_mat_raw[, meta_sub$SampleLabel, drop = FALSE]\ncount_mat <- round(count_mat)\nstorage.mode(count_mat) <- "integer"\n\nmeta_sub$genotype <- factor(meta_sub$genotype, levels = c("N2", "RB2060"))\nmeta_sub$treatment <- factor(meta_sub$treatment, levels = c("untreated", "treated"))\nmeta_sub$experiment <- factor(meta_sub$experiment)\n\ngene_symbols <- rownames(count_mat)\n\ncat("\\nFigure 5 design:\\n")\nprint(with(meta_sub, table(genotype, treatment, experiment)))\ncat("Samples:", nrow(meta_sub), "\\n")\ncat("Genes:", nrow(count_mat), "\\n")\n\ncat(\n  "\\nRunning model:\\n",\n  "~ genotype * treatment + genotype:experiment + treatment:experiment\\n"\n)\n\ndds <- DESeqDataSetFromMatrix(\n  countData = count_mat,\n  colData   = meta_sub,\n  design    = ~ genotype * treatment +\n              genotype:experiment +\n              treatment:experiment\n)\n\ndds <- DESeq(dds)\n\ncat("\\nCoefficient names:\\n")\nprint(resultsNames(dds))\n\nhave_apeglm <- requireNamespace("apeglm", quietly = TRUE)\nhave_ashr   <- requireNamespace("ashr", quietly = TRUE)\n\nif (!have_apeglm || !have_ashr) {\n  stop(\n    "Figure 5 requires both apeglm and ashr.\\n",\n    "Install with:\\n",\n    "  BiocManager::install(c(\'apeglm\',\'ashr\'))"\n  )\n}\n\nshrink_result <- function(dds, coef_name = NULL, contrast_list = NULL, type) {\n  if (!is.null(coef_name)) {\n    return(as.data.frame(lfcShrink(dds, coef = coef_name, type = type)))\n  }\n  as.data.frame(lfcShrink(dds, contrast = contrast_list, type = type))\n}\n\nres_treat_N2 <- shrink_result(\n  dds,\n  coef_name = "treatment_treated_vs_untreated",\n  type = "apeglm"\n)\nres_treat_N2$gene_symbol <- gene_symbols\n\nres_treat_RB <- shrink_result(\n  dds,\n  contrast_list = list(\n    c(\n      "treatment_treated_vs_untreated",\n      "genotypeRB2060.treatmenttreated"\n    )\n  ),\n  type = "ashr"\n)\nres_treat_RB$gene_symbol <- gene_symbols\n\nres_geno_unt <- shrink_result(\n  dds,\n  coef_name = "genotype_RB2060_vs_N2",\n  type = "apeglm"\n)\nres_geno_unt$gene_symbol <- gene_symbols\n\nres_geno_trt <- shrink_result(\n  dds,\n  contrast_list = list(\n    c(\n      "genotype_RB2060_vs_N2",\n      "genotypeRB2060.treatmenttreated"\n    )\n  ),\n  type = "ashr"\n)\nres_geno_trt$gene_symbol <- gene_symbols\n\nint_name <- grep(\n  "^genotypeRB2060\\\\.treatmenttreated$",\n  resultsNames(dds),\n  value = TRUE\n)\n\nif (length(int_name) != 1L) {\n  stop(\n    "Could not uniquely identify genotypeRB2060.treatmenttreated. Candidates: ",\n    paste(\n      grep(\n        "genotype.*RB2060.*treatment.*treated",\n        resultsNames(dds),\n        value = TRUE\n      ),\n      collapse = ", "\n    )\n  )\n}\n\nres_interaction <- shrink_result(\n  dds,\n  coef_name = int_name,\n  type = "apeglm"\n)\nres_interaction$gene_symbol <- gene_symbols\n\nwrite.csv(res_treat_N2, "DESeq2_treatment_in_N2.csv", row.names = TRUE)\nwrite.csv(res_treat_RB, "DESeq2_treatment_in_RB2060.csv", row.names = TRUE)\nwrite.csv(res_geno_unt, "DESeq2_genotype_untreated.csv", row.names = TRUE)\nwrite.csv(res_geno_trt, "DESeq2_genotype_treated.csv", row.names = TRUE)\nwrite.csv(res_interaction, "DESeq2_GxE_interaction.csv", row.names = TRUE)\n\nLFC_THRESH  <- 1.0\nPADJ_THRESH <- 0.05\n\nclassify <- data.frame(\n  gene_symbol   = gene_symbols,\n  lfc_treat_N2  = res_treat_N2$log2FoldChange,\n  padj_treat_N2 = res_treat_N2$padj,\n  lfc_treat_RB  = res_treat_RB$log2FoldChange,\n  padj_treat_RB = res_treat_RB$padj,\n  lfc_geno_unt  = res_geno_unt$log2FoldChange,\n  padj_geno_unt = res_geno_unt$padj,\n  lfc_geno_trt  = res_geno_trt$log2FoldChange,\n  padj_geno_trt = res_geno_trt$padj,\n  lfc_gxe       = res_interaction$log2FoldChange,\n  padj_gxe      = res_interaction$padj,\n  stringsAsFactors = FALSE\n)\n\nsig_tn2 <- !is.na(classify$padj_treat_N2) &\n  classify$padj_treat_N2 < PADJ_THRESH &\n  abs(classify$lfc_treat_N2) > LFC_THRESH\n\nsig_trb <- !is.na(classify$padj_treat_RB) &\n  classify$padj_treat_RB < PADJ_THRESH &\n  abs(classify$lfc_treat_RB) > LFC_THRESH\n\nsig_gut <- !is.na(classify$padj_geno_unt) &\n  classify$padj_geno_unt < PADJ_THRESH &\n  abs(classify$lfc_geno_unt) > LFC_THRESH\n\nsig_gtr <- !is.na(classify$padj_geno_trt) &\n  classify$padj_geno_trt < PADJ_THRESH &\n  abs(classify$lfc_geno_trt) > LFC_THRESH\n\nsig_gxe <- !is.na(classify$padj_gxe) &\n  classify$padj_gxe < PADJ_THRESH\n\nsig_treat <- sig_tn2 | sig_trb\nsig_geno  <- sig_gut | sig_gtr\n\nclassify$effect_class <- "NS"\nclassify$effect_class[sig_treat & !sig_geno] <- "Treatment"\nclassify$effect_class[!sig_treat & sig_geno] <- "Genotype"\nclassify$effect_class[sig_treat & sig_geno & sig_gxe] <- "GxE"\nclassify$effect_class[sig_treat & sig_geno & !sig_gxe] <- "Additive"\n\nwrite.csv(\n  classify,\n  "DESeq2_effect_classification.csv",\n  row.names = FALSE\n)\n\ncat("\\nEffect classification:\\n")\nprint(table(classify$effect_class))\n\ngxe <- classify[classify$effect_class == "GxE", , drop = FALSE]\n\nn_gxe <- nrow(gxe)\nn_sub <- sum(gxe$lfc_treat_RB <= gxe$lfc_treat_N2, na.rm = TRUE)\nn_supra <- sum(gxe$lfc_treat_RB > gxe$lfc_treat_N2, na.rm = TRUE)\n\ncat(\n  sprintf(\n    "\\nFigure 5 validation: %d GxE = %d sub + %d supra\\n",\n    n_gxe,\n    n_sub,\n    n_supra\n  )\n)\n\nif (n_gxe != 90L || n_sub != 68L || n_supra != 22L) {\n  warning(\n    "Figure 5 validation differs from recovered thesis result ",\n    "(90 GxE; 68 sub; 22 supra)."\n  )\n}\n'

def rebuild_effect_table() -> None:
    with tempfile.TemporaryDirectory(prefix="figure5_deseq2_") as tmp:
        tmpdir = Path(tmp)

        for src, name in (
            (RNA_FILE, "RNAseq_to_TF_Targets.csv"),
            (META_FILE, "Sample_Metadata_Table.csv"),
        ):
            dst = tmpdir / name
            if src.suffix == ".gz":
                # The R snippet expects a plain CSV; decompress on the way in.
                import gzip
                with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
                continue
            try:
                os.symlink(src.resolve(), dst)
            except OSError:
                shutil.copy2(src, dst)

        r_file = tmpdir / "Figure5_effect_classification.R"
        r_file.write_text(R_EFFECT_SCRIPT)

        print("[analysis] Rebuilding Figure 5 DESeq2 effect table...")
        result = subprocess.run(["Rscript", str(r_file)], cwd=tmpdir)

        if result.returncode != 0:
            sys.exit(
                f"[error] Figure 5 DESeq2 analysis failed "
                f"(exit code {result.returncode})."
            )

        expected = (
            "DESeq2_treatment_in_N2.csv",
            "DESeq2_treatment_in_RB2060.csv",
            "DESeq2_genotype_untreated.csv",
            "DESeq2_genotype_treated.csv",
            "DESeq2_GxE_interaction.csv",
            "DESeq2_effect_classification.csv",
        )

        for filename in expected:
            src = tmpdir / filename
            if not src.exists():
                sys.exit(f"[error] Expected R output missing: {filename}")
            shutil.copy2(src, DESEQ_DIR / filename)

    print(f"[ok] Rebuilt effect table: {EFFECT_FILE}")

PADJ_CUT = 0.05

HARDCODED_NAMES = {
    "B0284.2": "pals-27", "C10C5.2": "fbxc-58",
    "C29F9.4": "pals-24", "C31A11.7": "oac-7",
    "C45G7.3": "ilys-3", "F08F3.7": "cyp-14a5",
    "F28D1.3": "thn-1", "F35C5.10": "nspb-11",
    "F37B1.8": "gst-19", "K10G6.4": "cank-26",
    "R10D12.17": "srw-145", "T06E6.5": "fbxa-135",
    "T08G5.10": "mtl-2", "T21E8.2": "pgp-7",
    "Y39G10AR.6": "ugt-31", "Y47D3A.2": "fbxa-128",
    "Y47D3B.10": "dpy-18", "Y49C4A.8": "ugt-29",
}

PANEL_C_TFS = [
    "blmp-1","cebp-1","ceh-60","daf-16","elt-2","fos-1",
    "lin-35","nhr-28","nhr-77","pha-4","pqm-1","sma-9","snpc-4","zip-2",
]

COLOR_SUB = "#2471A3"
COLOR_SUPRA = "#C0392B"
HEATMAP_VMAX = 4.0
SCATTER_LIM = 9.0
N_SUB_SELECT = 30
N_SUPRA_SELECT = 8
MIN_SUPRA_LFC = 1.0

COL_GENO_UNT = "lfc_geno_unt"
COL_TREAT_N2 = "lfc_treat_N2"
COL_GENO_TRT = "lfc_geno_trt"
COL_GXE = "lfc_gxe"
COL_TREAT_RB = "lfc_treat_RB"

CAT_COLORS = {
    "Treatment-only": "#F8766D",
    "Genotype-only": "#00BA38",
    "Conditional interaction": "#619CFF",
}

FS_PANEL = 9
FS_TITLE = 7
FS_LABEL = 7
FS_TICK = 6
FS_ANNOT = 6

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS_TICK,
    "axes.titlesize": FS_TITLE,
    "axes.titleweight": "bold",
    "axes.labelsize": FS_LABEL,
    "xtick.labelsize": FS_TICK,
    "ytick.labelsize": FS_TICK,
    "legend.fontsize": FS_ANNOT,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.0,
    "ytick.major.size": 2.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.default": "regular",
})

MM = 1 / 25.4

def load_effects() -> pd.DataFrame:
    genome = pd.read_csv(EFFECT_FILE)
    numeric_cols = [
        COL_TREAT_N2, COL_TREAT_RB, COL_GENO_UNT, COL_GENO_TRT, COL_GXE,
        "padj_treat_N2", "padj_treat_RB", "padj_geno_unt",
        "padj_geno_trt", "padj_gxe",
    ]
    for col in numeric_cols:
        genome[col] = pd.to_numeric(genome[col], errors="coerce")
    return genome

def build_name_lookup(genome: pd.DataFrame) -> dict[str, str]:
    lookup = dict(HARDCODED_NAMES)
    merged = pd.read_csv(RNA_FILE, low_memory=False)

    if {"Gene", "Name.Target"}.issubset(merged.columns):
        for _, row in merged[["Gene", "Name.Target"]].drop_duplicates().iterrows():
            gene = str(row["Gene"]).strip()
            name = str(row["Name.Target"]).strip()
            if gene and name not in ("-", "", "nan"):
                lookup[gene] = name

    return lookup

def select_panel_b_genes(genome: pd.DataFrame) -> pd.DataFrame:
    gxe = genome[genome["effect_class"] == "GxE"].copy()
    gxe = gxe.dropna(subset=[COL_GXE])

    sub = gxe[gxe[COL_GXE] < 0].copy()
    sub["_score"] = -sub[COL_GXE] + sub[COL_GENO_UNT].fillna(0).clip(lower=0)
    sub_sel = sub.nlargest(N_SUB_SELECT, "_score")

    supra = gxe[gxe[COL_GXE] >= MIN_SUPRA_LFC].copy()
    supra_sel = supra.nlargest(N_SUPRA_SELECT, COL_GXE)

    selected = (
        pd.concat([sub_sel, supra_sel])
        .drop_duplicates("gene_symbol")
        .drop(columns=["_score"], errors="ignore")
    )

    print(
        f"[ok] Panel B: {len(sub_sel)} sub-additive + "
        f"{len(supra_sel)} supra-additive = {len(selected)} genes"
    )
    return selected

def build_panel_c_data(genome: pd.DataFrame) -> pd.DataFrame:
    merged = pd.read_csv(RNA_FILE, low_memory=False)
    tf_map = merged[["Gene", "Name.TF"]].dropna().drop_duplicates()
    tf_map["Gene"] = tf_map["Gene"].astype(str).str.strip()
    tf_map["Name.TF"] = tf_map["Name.TF"].astype(str).str.lower().str.strip()

    def sig_set(padj_col: str) -> set[str]:
        return set(
            genome.loc[
                genome[padj_col].notna() & (genome[padj_col] < PADJ_CUT),
                "gene_symbol",
            ].astype(str)
        )

    de_treat_n2 = sig_set("padj_treat_N2")
    de_treat_rb = sig_set("padj_treat_RB")
    de_geno_unt = sig_set("padj_geno_unt")
    de_geno_trt = sig_set("padj_geno_trt")

    rows = []
    for tf_name in PANEL_C_TFS:
        targets = set(
            tf_map.loc[tf_map["Name.TF"] == tf_name, "Gene"]
        )
        total = len(targets)

        if total == 0:
            rows.append({
                "TF": tf_name.upper(),
                "Treatment-only": 0.0,
                "Genotype-only": 0.0,
                "Conditional interaction": 0.0,
            })
            continue

        treat_union = (de_treat_n2 | de_treat_rb) & targets
        geno_union = (de_geno_unt | de_geno_trt) & targets

        rows.append({
            "TF": tf_name.upper(),
            "Treatment-only": len(treat_union - geno_union) / total,
            "Genotype-only": len(geno_union - treat_union) / total,
            "Conditional interaction": len(treat_union & geno_union) / total,
        })

    return pd.DataFrame(rows)

def panel_a(ax, genome: pd.DataFrame, panel_b_genes: set[str]) -> None:
    df = genome[
        genome[COL_TREAT_N2].notna()
        & genome[COL_TREAT_RB].notna()
        & genome["effect_class"].isin(["GxE", "Treatment", "Additive", "Genotype"])
    ].copy()

    gxe = df[df["effect_class"] == "GxE"].copy()
    lim = SCATTER_LIM

    ax.plot([-lim, lim], [-lim, lim], color="#555", lw=0.8,
            label="additive (y=x)", zorder=1)
    for offset in (1, -1):
        ax.plot([-lim, lim], [-lim + offset, lim + offset],
                color="#bbb", lw=0.35, ls=":", zorder=1)

    other_response = df[df["effect_class"] != "GxE"]
    ax.scatter(other_response[COL_TREAT_N2].clip(-lim, lim),
               other_response[COL_TREAT_RB].clip(-lim, lim),
               s=1.5, c="#BBBBBB", alpha=0.30, edgecolors="none",
               label="other responsive", zorder=2)

    other_gxe = gxe[~gxe["gene_symbol"].isin(panel_b_genes)]
    ax.scatter(other_gxe[COL_TREAT_N2].clip(-lim, lim),
               other_gxe[COL_TREAT_RB].clip(-lim, lim),
               s=8, c="#E07B3D", alpha=0.40, edgecolors="none",
               label=f"G×E not in b (n={len(other_gxe)})", zorder=3)

    selected = gxe[gxe["gene_symbol"].isin(panel_b_genes)].copy()
    selected["above"] = selected[COL_TREAT_RB] > selected[COL_TREAT_N2]
    selected_sub = selected[~selected["above"]]
    selected_supra = selected[selected["above"]]

    ax.scatter(selected_sub[COL_TREAT_N2].clip(-lim, lim),
               selected_sub[COL_TREAT_RB].clip(-lim, lim),
               s=22, c=COLOR_SUB, edgecolors="white", linewidths=0.3,
               label=f"sub-additive in b (n={len(selected_sub)})", zorder=5)

    ax.scatter(selected_supra[COL_TREAT_N2].clip(-lim, lim),
               selected_supra[COL_TREAT_RB].clip(-lim, lim),
               s=22, c=COLOR_SUPRA, edgecolors="white", linewidths=0.3,
               label=f"supra-additive in b (n={len(selected_supra)})", zorder=5)

    n_gxe = len(gxe)
    n_sub = int((gxe[COL_TREAT_RB] <= gxe[COL_TREAT_N2]).sum())
    n_supra = int((gxe[COL_TREAT_RB] > gxe[COL_TREAT_N2]).sum())
    k = min(n_sub, n_supra)
    binom_p = min(
        1.0,
        2.0 * sum(comb(n_gxe, i) for i in range(k + 1)) / (2 ** n_gxe),
    )

    ax.text(
        -lim + 0.3, lim - 0.3,
        f"{n_gxe} G×E genes:\n"
        f"{n_sub} sub-additive, {n_supra} supra-additive\n"
        f"binomial p = {binom_p:.2g}",
        ha="left", va="top", fontsize=FS_ANNOT,
        bbox=dict(boxstyle="round,pad=0.2", fc="white",
                  ec="#888", lw=0.5, alpha=0.9),
    )

    ax.text(lim - 0.3, lim - 0.3, "supra-additive",
            ha="right", va="top", fontsize=FS_ANNOT,
            style="italic", color=COLOR_SUPRA,
            bbox=dict(boxstyle="round,pad=0.18", fc="white",
                      ec=COLOR_SUPRA, lw=0.5, alpha=0.85))

    ax.text(lim - 0.3, -lim + 0.3, "sub-additive",
            ha="right", va="bottom", fontsize=FS_ANNOT,
            style="italic", color=COLOR_SUB,
            bbox=dict(boxstyle="round,pad=0.18", fc="white",
                      ec=COLOR_SUB, lw=0.5, alpha=0.85))

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xticks([-8, -4, 0, 4, 8])
    ax.set_yticks([-8, -4, 0, 4, 8])
    ax.set_xlabel("CoCl$_2$ response in N2 (log$_2$FC)")
    ax.set_ylabel("CoCl$_2$ response in argk-2$^{-/-}$ (log$_2$FC)")
    ax.set_title("Non-additive CoCl$_2$ response in argk-2$^{-/-}$", pad=3)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.legend(loc="lower right", fontsize=FS_ANNOT, frameon=True,
              framealpha=0.90, edgecolor="#aaa",
              handlelength=1.0, handletextpad=0.35,
              borderpad=0.35, labelspacing=0.3)

    ax.text(-0.20, 1.04, "a", transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold")

def panel_b(ax, ax_cb, selected: pd.DataFrame,
            name_lookup: dict[str, str]) -> None:
    sub = selected[selected[COL_GXE] < 0].sort_values(COL_GXE, ascending=True)
    supra = selected[selected[COL_GXE] >= MIN_SUPRA_LFC].sort_values(
        COL_GXE, ascending=False
    )
    effectors = pd.concat([sub, supra]).copy()
    n_sub = len(sub)

    effectors["display_name"] = (
        effectors["gene_symbol"].astype(str)
        .map(lambda gene: name_lookup.get(gene, gene))
    )

    matrix = np.clip(
        effectors[[COL_GENO_UNT, COL_TREAT_N2, COL_GENO_TRT, COL_GXE]]
        .fillna(0).values,
        -HEATMAP_VMAX,
        HEATMAP_VMAX,
    )

    n_genes, n_cols = matrix.shape

    ax.imshow(matrix, cmap=plt.cm.RdBu_r, aspect="auto",
              vmin=-HEATMAP_VMAX, vmax=HEATMAP_VMAX,
              extent=[0, n_cols, 0, n_genes],
              origin="upper", interpolation="nearest")

    ax.plot([3, 3], [0, n_genes], color="black", lw=1.0, zorder=4)

    if 0 < n_sub < n_genes:
        ax.plot([0, n_cols], [n_sub, n_sub], color="#333",
                lw=0.8, ls="--", zorder=4)

        # origin="upper": first block is displayed at the top.
        ax.text(n_cols + 0.1, n_sub / 2, "sub-additive",
                va="center", ha="left", fontsize=FS_ANNOT,
                style="italic", color=COLOR_SUB, clip_on=False)

        ax.text(n_cols + 0.1, n_sub + (n_genes - n_sub) / 2,
                "supra-additive",
                va="center", ha="left", fontsize=FS_ANNOT,
                style="italic", color=COLOR_SUPRA, clip_on=False)

    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(
        ["Genotype", "Treatment", "Genotype\n(CoCl$_2$)", "G×E\ninteraction"],
        fontsize=FS_TICK, ha="center", linespacing=1.1,
    )
    ax.xaxis.set_ticks_position("top")

    ax.set_yticks(np.arange(n_genes) + 0.5)
    ax.set_yticklabels(effectors["display_name"].tolist(),
                       fontsize=FS_TICK, style="italic")
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", length=0, pad=2)

    for spine in ("top", "right", "left", "bottom"):
        ax.spines[spine].set_visible(False)

    ax.set_title(
        "G×E effector expression\n(sorted by interaction magnitude)",
        pad=26,
    )
    ax.text(-0.42, 1.04, "b", transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold")

    mapper = plt.cm.ScalarMappable(
        cmap=plt.cm.RdBu_r,
        norm=mpl.colors.Normalize(-HEATMAP_VMAX, HEATMAP_VMAX),
    )
    colorbar = plt.colorbar(
        mapper, cax=ax_cb, orientation="horizontal",
        ticks=[-HEATMAP_VMAX, 0, HEATMAP_VMAX],
    )
    colorbar.set_label("log$_2$FC", fontsize=FS_TICK, labelpad=1)
    colorbar.ax.tick_params(labelsize=FS_TICK, length=1.5, width=0.5)
    colorbar.outline.set_linewidth(0.4)

def panel_c(ax, proportions: pd.DataFrame) -> None:
    categories = list(CAT_COLORS.keys())
    x = np.arange(len(proportions))
    bottom = np.zeros(len(proportions))

    for category in categories:
        values = proportions[category].values
        ax.bar(x, values, bottom=bottom, color=CAT_COLORS[category],
               label=category, width=0.72, linewidth=0, zorder=3)
        bottom += values

    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#dddddd", linewidth=0.5, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(proportions["TF"].tolist(),
                       rotation=45, ha="right", fontsize=FS_TICK)
    ax.set_ylabel("Proportion of TF targets")
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
    ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1.00"],
                       fontsize=FS_TICK)
    ax.set_title("Effect category of TF targets", pad=3)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.legend(loc="upper right", fontsize=FS_ANNOT,
              frameon=True, framealpha=0.92,
              edgecolor="#cccccc", facecolor="white",
              title="Effect Category", title_fontsize=FS_ANNOT)

    ax.text(-0.20, 1.04, "c", transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold")

def main() -> None:
    rebuild_effect_table()

    genome = load_effects()
    name_lookup = build_name_lookup(genome)
    selected = select_panel_b_genes(genome)
    panel_b_genes = set(selected["gene_symbol"].dropna().astype(str))
    panel_c_df = build_panel_c_data(genome)

    selected.assign(
        display_name=selected["gene_symbol"].astype(str)
        .map(lambda gene: name_lookup.get(gene, gene))
    ).to_csv(
        FIG_DIR / "Figure5_panelB_selected_genes.csv",
        index=False,
    )

    genome.to_csv(
        FIG_DIR / "Figure5_effect_table_used.csv",
        index=False,
    )

    n_genes = len(selected)
    heatmap_data_h_mm = n_genes * 3.4
    fig_h_mm = min(max(heatmap_data_h_mm + 24, 120), 230)

    fig_w = 183 * MM
    fig_h = fig_h_mm * MM
    fig = plt.figure(figsize=(fig_w, fig_h))

    left_col_w_mm = 183 * 0.38
    scatter_h_frac = left_col_w_mm / fig_h_mm
    bar_h_frac = 1.0 - scatter_h_frac - 0.04
    cb_h_frac = 6 / fig_h_mm

    outer = fig.add_gridspec(
        1, 2,
        width_ratios=[0.46, 0.54],
        wspace=0.26,
        left=0.10, right=0.88,
        top=0.95, bottom=0.06,
    )

    left = outer[0, 0].subgridspec(
        3, 1,
        height_ratios=[scatter_h_frac, 0.03, bar_h_frac],
        hspace=0.0,
    )
    ax_a = fig.add_subplot(left[0, 0])
    ax_c = fig.add_subplot(left[2, 0])

    right = outer[0, 1].subgridspec(
        2, 1,
        height_ratios=[1.0 - cb_h_frac, cb_h_frac],
        hspace=0.04,
    )
    ax_b = fig.add_subplot(right[0, 0])
    ax_cb = fig.add_subplot(right[1, 0])

    panel_a(ax_a, genome, panel_b_genes)
    panel_b(ax_b, ax_cb, selected, name_lookup)
    panel_c(ax_c, panel_c_df)

    plt.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    plt.savefig(OUT_PNG, dpi=600, bbox_inches="tight")
    plt.close(fig)

    gxe = genome[genome["effect_class"] == "GxE"]
    n_sub = int((gxe[COL_TREAT_RB] <= gxe[COL_TREAT_N2]).sum())
    n_supra = int((gxe[COL_TREAT_RB] > gxe[COL_TREAT_N2]).sum())

    print("\n========================================")
    print("FIGURE 5 COMPLETE")
    print("========================================")
    print(f"G×E genes: {len(gxe)}")
    print(f"sub-additive: {n_sub}")
    print(f"supra-additive: {n_supra}")
    print(f"PDF: {OUT_PDF}")
    print(f"PNG: {OUT_PNG}")

if __name__ == "__main__":
    main()
