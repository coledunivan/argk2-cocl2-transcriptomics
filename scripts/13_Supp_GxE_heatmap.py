#!/usr/bin/env python3
"""
Supplementary Figure 3 — full G×E vertical heatmap, no TF-module grouping.

Renders all 90 G×E genes as rows against the four DESeq2 contrasts as columns:

    1. Genotype        argk-2−/− untreated vs N2 untreated
    2. Treatment       N2 +CoCl2 vs N2 untreated
    3. Stress genotype argk-2−/− +CoCl2 vs N2 +CoCl2
    4. G×E             non-additive interaction term

Genes are sorted by the G×E interaction value, so the sub-additive and
supra-additive blocks separate visually without needing module labels.
Every G×E gene is retained, including those with no TFLink regulator.

INPUTS
    outputs/DESeq2/DESeq2_effect_classification.csv   (from 01_DESeq2_and_TF_enrichment.R)
    data/reference/RNAseq_to_TF_Targets.csv[.gz]      (common gene-name lookup only)

OUTPUTS
    outputs/supplementary/Figure_S3_GxE_heatmap.pdf
    outputs/supplementary/Figure_S3_GxE_heatmap.png

Usage:
    python3 scripts/14_Supp_GxE_heatmap.py
"""

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- Repository-anchored paths ---------------------------------------------
# Resolve all inputs and outputs from the repo root so the script can be
# launched from anywhere. Expected layout: <repo>/scripts/<script>.py
SCRIPT_PATH = Path(__file__).resolve()
ROOT = SCRIPT_PATH.parent.parent if SCRIPT_PATH.parent.name == "scripts" else Path.cwd()

DESEQ_FILE = ROOT / "outputs" / "DESeq2" / "DESeq2_effect_classification.csv"
OUT_DIR = ROOT / "outputs" / "supplementary"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF = OUT_DIR / "Figure_S3_GxE_heatmap.pdf"
OUT_PNG = OUT_DIR / "Figure_S3_GxE_heatmap.png"


def resolve_tflink():
    """Prefer the plain CSV; fall back to the gzipped copy shipped on GitHub."""
    for name in ("RNAseq_to_TF_Targets.csv", "RNAseq_to_TF_Targets.csv.gz"):
        p = ROOT / "data" / "reference" / name
        if p.exists():
            return p
    return None


TFLINK_FILE = resolve_tflink()

# Optional shared contrast-validation helper. Absent from this repository, so
# the import is guarded; the figure is identical either way.
try:
    from _contrasts import validate_contrast_columns, provenance_footer
except ImportError:
    validate_contrast_columns = None
    provenance_footer = None

HEATMAP_VMAX = 4.0

# Actual column names in the DESeq2 classification table.
# Order here controls left-to-right order in the heatmap.
HEATMAP_COLS = [
    "lfc_geno_unt",   # argk-2 untreated vs N2 untreated
    "lfc_treat_N2",   # N2 +CoCl2 vs N2 untreated
    "lfc_geno_trt",   # argk-2 +CoCl2 vs N2 +CoCl2
    "lfc_gxe",        # non-additive interaction term
]

HEATMAP_LABELS = ["Genotype", "Treatment", "Stress\ngenotype", "G×E"]

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 0.9,
    "pdf.fonttype": 42,
})


def load_gxe_genes():
    """Load the G×E gene set from the DESeq2 effect classification table."""
    if not DESEQ_FILE.exists():
        sys.exit(
            f"[fatal] {DESEQ_FILE} not found.\n"
            "        Run 01_DESeq2_and_TF_enrichment.R first."
        )

    print(f"[ok] reading G×E classification from {DESEQ_FILE.relative_to(ROOT)}")
    genome = pd.read_csv(DESEQ_FILE)

    for col in HEATMAP_COLS + ["padj_gxe", "lfc_treat_RB"]:
        if col in genome.columns:
            genome[col] = pd.to_numeric(genome[col], errors="coerce")

    missing = [c for c in HEATMAP_COLS if c not in genome.columns]
    if missing:
        sys.exit(f"[fatal] Missing required contrast columns: {missing}")

    if validate_contrast_columns is not None:
        validate_contrast_columns(genome)

    gxe = genome[genome["effect_class"] == "GxE"].copy()
    print(f"     {len(gxe)} G×E genes from DESeq2")
    return gxe, genome


def load_tflink_for_names():
    """Load the TFLink table for common-name lookup only."""
    if TFLINK_FILE is None:
        print("[warn] TFLink table not found. Using gene IDs as labels.")
        return None
    print(f"[ok] reading gene-name lookup from {TFLINK_FILE.relative_to(ROOT)}")
    # pandas decompresses .gz transparently from the file extension.
    return pd.read_csv(TFLINK_FILE, low_memory=False)


def build_name_lookup(tf_df):
    """One pass over TFLink -> {gene_id: common_name}. Avoids a scan per gene."""
    if tf_df is None or "Gene" not in tf_df.columns or "Name.Target" not in tf_df.columns:
        return {}
    sub = tf_df[["Gene", "Name.Target"]].dropna()
    sub = sub[sub["Name.Target"] != "-"]
    sub = sub.drop_duplicates(subset="Gene", keep="first")
    return {
        g: str(n).split(";")[0].split(",")[0].strip().lower()
        for g, n in zip(sub["Gene"], sub["Name.Target"])
    }


def main():
    gxe, _ = load_gxe_genes()
    tf_raw = load_tflink_for_names()
    name_lookup = build_name_lookup(tf_raw)

    gene_col = "gene_symbol" if "gene_symbol" in gxe.columns else "Gene"
    if gene_col not in gxe.columns:
        sys.exit("[fatal] No gene ID column: expected 'gene_symbol' or 'Gene'.")

    rows = []
    for _, row in gxe.iterrows():
        gene_id = row[gene_col]
        rows.append({
            "gene_id": gene_id,
            "label": name_lookup.get(gene_id, gene_id),
            **{k: row.get(k, np.nan) for k in HEATMAP_COLS},
        })
    df = pd.DataFrame(rows)

    # No TF/module grouping. Sorting by the interaction term preserves the
    # sub-additive -> supra-additive structure without module labels.
    df = df.sort_values("lfc_gxe", ascending=True).reset_index(drop=True)

    mat = np.clip(df[HEATMAP_COLS].fillna(0).values, -HEATMAP_VMAX, HEATMAP_VMAX)
    n_genes = len(df)
    print(f"\nHeatmap: {n_genes} G×E genes, no TF-module grouping")

    fig_h = max(9, n_genes * 0.18)
    fig, ax = plt.subplots(figsize=(5.8, fig_h))

    im = ax.imshow(mat, cmap=plt.cm.RdBu_r, aspect="auto",
                   vmin=-HEATMAP_VMAX, vmax=HEATMAP_VMAX, interpolation="nearest")

    ax.set_xticks(np.arange(len(HEATMAP_COLS)))
    ax.set_xticklabels(HEATMAP_LABELS, fontsize=10)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=6)

    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels(df["label"], fontsize=7)
    ax.tick_params(axis="y", length=0, pad=2)

    ax.set_xticks(np.arange(-0.5, len(HEATMAP_COLS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_genes, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)

    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    ax.set_title(
        f"Directional response of all {n_genes} G×E genes\n"
        "Genotype = argk-2−/− untreated vs N2 untreated; "
        "Treatment = N2 +CoCl₂ vs N2 untreated; "
        "Stress genotype = argk-2−/− +CoCl₂ vs N2 +CoCl₂",
        fontsize=10.5, fontweight="bold", pad=36,
    )

    cax = fig.add_axes([0.86, 0.25, 0.035, 0.5])
    cb = fig.colorbar(im, cax=cax, orientation="vertical",
                      ticks=[-HEATMAP_VMAX, 0, HEATMAP_VMAX])
    cb.set_label("log$_2$FC", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    fig.text(0.12, 0.02,
             "G×E = genotype-by-treatment interaction term; "
             "genes sorted by G×E effect size.",
             ha="left", va="bottom", fontsize=7, color="0.35")

    if provenance_footer is not None:
        provenance_footer(fig, contrast_keys=HEATMAP_COLS,
                          script_name=SCRIPT_PATH.name)

    plt.savefig(OUT_PDF, dpi=400, bbox_inches="tight")
    plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"\n[done] {OUT_PDF.relative_to(ROOT)}  +  {OUT_PNG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
