#!/usr/bin/env python3
"""Rebuild the main figure from the already-saved Pfaffl-corrected CSVs.

Faster than rerunning apply_pfaffl.py because it doesn't need to re-do the
Cq loading / QC / DeltaCq / aggregation steps. Use this whenever the
analysis math hasn't changed but figure styling has (labels, colors, etc.).
"""
import sys
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "qpcr"))
from qpcr_analysis import (
    make_figure, MIN_RUNS_FOR_FIG, GENES_EXCLUDE, QPCR_OUT, SUPP_FIG,
)

# Read the Pfaffl-corrected tables produced by apply_pfaffl.py
agg = pd.read_csv(QPCR_OUT / "qPCR_combined_deltadeltaCq_pfaffl.csv")
per_run = pd.read_csv(QPCR_OUT / "qPCR_per_run_ddcq_pfaffl.csv")

# Apply the same figure-subset filtering used in main()
agg_fig = agg[agg["n_runs"] >= MIN_RUNS_FOR_FIG].copy()
genes_with_key = set(agg_fig[agg_fig["contrast"] == "geno_trt"]["target"])
agg_fig = agg_fig[agg_fig["target"].isin(genes_with_key) &
                  ~agg_fig["target"].isin(GENES_EXCLUDE)]
ddcq_fig = per_run[per_run["target"].isin(genes_with_key) &
                   ~per_run["target"].isin(GENES_EXCLUDE)]

print(f"Figure subset: {len(agg_fig)} rows, {len(genes_with_key)} genes")
make_figure(agg_fig, ddcq_fig,
            SUPP_FIG / "qPCR_figure_panel.pdf",
            SUPP_FIG / "qPCR_figure_panel.png")
