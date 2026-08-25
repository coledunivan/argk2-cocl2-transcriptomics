#!/usr/bin/env python3
"""
add_new_run.py - Add ColeD_41726 qPCR run to the combined analysis.

Strategy
The existing CSVs (qPCR_per_run_ddcq.csv and qPCR_combined_deltadeltaCq.csv)
summarize the 7 previously processed plates. 
  1. Load the existing per-run table (old runs; 7 plates).
  2. Process ONLY the new run (ColeD_41726) through the same pipeline
     (outlier-flag -> DeltaCq vs tba-1 -> DeltaDeltaCq for the 3 contrasts).
  3. Concatenate old + new per-run rows.
  4. Re-run aggregate_across_runs() over the union so the new run
     contributes to the mean, SEM, and one-sample t-test of every
     (gene, contrast) row.
  5. Refresh qPCR_vs_RNAseq_comparison.csv by replacing its qPCR columns
     with the new aggregate (RNA-seq columns are untouched because the
     RNA-seq DESeq2 results have not changed).
  6. Rebuild qPCR_figure_panel.pdf/png via make_figure().

"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qpcr_analysis import (
    load_all_cq, flag_outliers, compute_delta_cq, compute_ddcq_contrasts,
    aggregate_across_runs, bh_adjust, make_figure,
    TARGET_RENAMES, REF_GENE, MIN_RUNS_FOR_FIG, GENES_EXCLUDE,
)

from qpcr_analysis import QPCR_OUT, QPCR_ROOT, SUPP_FIG  # repo-anchored paths

WORK = QPCR_ROOT
NEW_RUN_DIR = WORK / "ColeD_41726"
PER_RUN_CSV = QPCR_OUT / "qPCR_per_run_ddcq.csv"
COMBINED_CSV = QPCR_OUT / "qPCR_combined_deltadeltaCq.csv"
COMPARISON_CSV = QPCR_OUT / "qPCR_vs_RNAseq_comparison.csv"

# 1. Process the NEW run only
print("=" * 72)
print("Adding new run:", NEW_RUN_DIR.name)
print("=" * 72)

# load_all_cq iterates root.iterdir(), so we pass a parent that has only
# this one subdirectory. Our WORK dir has exactly one run folder, so that's
# what we want.
cq_df = load_all_cq(WORK)
print(f"\n[load] wells from new run: {len(cq_df)}")
print(f"       runs found:           {sorted(cq_df['run'].unique())}")
print(f"       targets:              {sorted(cq_df['target'].unique())}")
print(f"       strains:              {sorted(cq_df['strain'].unique())}")
print(f"       treatments:           {sorted(cq_df['treatment'].unique())}")
print(f"       condition counts:")
print(cq_df.groupby(['strain','treatment','target'])
      .size().unstack(fill_value=0).to_string())

# 2. QC: flag technical-replicate outliers (SD>0.5 drop logic)
print("\n[QC] flagging technical-replicate outliers (SD>0.5)")
cq_df = flag_outliers(cq_df)
dropped = cq_df[~cq_df['keep']]
print(f"     wells dropped: {len(dropped)}")
if len(dropped):
    print(dropped[['run','target','strain','treatment','well','cq']]
          .to_string(index=False))

# 3. DeltaCq vs tba-1

print("\n[DeltaCq] normalizing to tba-1")
dcq_df = compute_delta_cq(cq_df, ref_gene=REF_GENE)
print(f"          rows: {len(dcq_df)}")

# 4. DeltaDeltaCq for 3 contrasts
print("\n[DeltaDeltaCq] computing 3 contrasts")
ddcq_new = compute_ddcq_contrasts(dcq_df)
# Apply the canonical rename (t22f3.3 -> t22f3.11). For this run the target
# was already named t22f3.11 on the plate, but we apply it unconditionally
# for safety.
ddcq_new["target"] = ddcq_new["target"].replace(TARGET_RENAMES)
ddcq_save = ddcq_new.drop(columns=["dcq_num_vals", "dcq_den_vals"])
print(ddcq_save[['run','target','contrast','ddcq','log2fc','n_num','n_den']]
      .to_string(index=False))

# 5. Load existing per-run CSV, concatenate, re-aggregate
print("\n[merge] concatenating with existing per-run table")
old_per_run = pd.read_csv(PER_RUN_CSV)
print(f"        old runs: {old_per_run['run'].nunique()} "
      f"({old_per_run['run'].nunique()} plates, {len(old_per_run)} rows)")

# Make sure we don't double-count if the new run was already merged
if NEW_RUN_DIR.name in set(old_per_run["run"]):
    print(f"        NEW RUN {NEW_RUN_DIR.name!r} already in old table - removing first")
    old_per_run = old_per_run[old_per_run["run"] != NEW_RUN_DIR.name]

merged_per_run = pd.concat(
    [old_per_run, ddcq_save[old_per_run.columns]],
    ignore_index=True,
)
print(f"        merged:   {merged_per_run['run'].nunique()} runs, "
      f"{len(merged_per_run)} rows")

# For aggregate_across_runs we need the original ddcq_df shape (with n_num,
# n_den columns). The old per-run CSV already has those, so we're good.
print("\n[aggregate] re-aggregating across all runs")
agg_new = aggregate_across_runs(merged_per_run)
agg_new["padj_BH"] = bh_adjust(agg_new["p_value"].values)
agg_new = agg_new.sort_values(["contrast", "target"]).reset_index(drop=True)

# 6. Save refreshed CSVs
merged_per_run.to_csv(PER_RUN_CSV, index=False)
print(f"\n[write] {PER_RUN_CSV}  ({len(merged_per_run)} rows)")

agg_new.to_csv(COMBINED_CSV, index=False)
print(f"[write] {COMBINED_CSV}  ({len(agg_new)} rows)")

# 7. Refresh qPCR_vs_RNAseq_comparison.csv
#    - keep rnaseq_* columns as-is (RNA-seq data hasn't changed)
#    - overwrite qpcr_* columns from the new aggregate
#    - recompute direction_agrees
print(f"\n[write] refreshing {COMPARISON_CSV}")
comp_old = pd.read_csv(COMPARISON_CSV)

# Build a (gene, contrast) -> new qPCR stats map
qpcr_map = agg_new.set_index(["target", "contrast"]).to_dict("index")

def _pick(gene, contrast, key):
    entry = qpcr_map.get((gene, contrast))
    return entry.get(key) if entry else np.nan

# Collapse the duplicate t22f3.11 / geno_unt row in the old comp CSV: the
# refreshed aggregate has exactly one row per (gene, contrast) so the
# duplicated row in the old comparison table is no longer meaningful.
comp_new = comp_old.drop_duplicates(subset=["gene", "contrast"]).copy()

comp_new["qpcr_log2fc"]  = [_pick(g, c, "log2fc_mean")
                            for g, c in zip(comp_new["gene"], comp_new["contrast"])]
comp_new["qpcr_sem"]     = [_pick(g, c, "log2fc_sem")
                            for g, c in zip(comp_new["gene"], comp_new["contrast"])]
comp_new["qpcr_n_runs"]  = [_pick(g, c, "n_runs")
                            for g, c in zip(comp_new["gene"], comp_new["contrast"])]
comp_new["qpcr_pvalue"]  = [_pick(g, c, "p_value")
                            for g, c in zip(comp_new["gene"], comp_new["contrast"])]
comp_new["qpcr_padj_BH"] = [_pick(g, c, "padj_BH")
                            for g, c in zip(comp_new["gene"], comp_new["contrast"])]

# Recompute direction_agrees
sign_q = np.sign(comp_new["qpcr_log2fc"])
sign_r = np.sign(comp_new["rnaseq_log2fc"])
both = comp_new["qpcr_log2fc"].notna() & comp_new["rnaseq_log2fc"].notna()
comp_new["direction_agrees"] = np.where(both, sign_q == sign_r, pd.NA)

comp_new.to_csv(COMPARISON_CSV, index=False)
print(f"        wrote {len(comp_new)} rows")

# 8. Report the effect of the new run on the figure-subset genes
print("\n" + "=" * 72)
print("EFFECT OF NEW RUN ON FIGURE-SUBSET GENES")
print("=" * 72)
FIG_GENES = ["ugt-31", "oac-7", "ugt-29", "t22f3.11"]

for gene in FIG_GENES:
    print(f"\n--- {gene} ---")
    for contrast in ["treat_N2", "geno_unt", "geno_trt"]:
        new_row = ddcq_save[(ddcq_save["target"] == gene) &
                            (ddcq_save["contrast"] == contrast)]
        agg_row = agg_new[(agg_new["target"] == gene) &
                          (agg_new["contrast"] == contrast)]
        if new_row.empty or agg_row.empty:
            continue
        this_run_fc = new_row["log2fc"].iloc[0]
        mean_fc = agg_row["log2fc_mean"].iloc[0]
        sem = agg_row["log2fc_sem"].iloc[0]
        n = int(agg_row["n_runs"].iloc[0])
        p = agg_row["p_value"].iloc[0]
        print(f"  {contrast:10s}  new-run log2FC = {this_run_fc:+.2f}   "
              f"combined: {mean_fc:+.2f} ± {sem:.2f}  n={n}  p={p:.3g}")

# 9. Build the figure
print("\n" + "=" * 72)
print("BUILDING FIGURE")
print("=" * 72)

# Reproduce the figure-filtering logic from qpcr_analysis.main()
agg_fig = agg_new[agg_new["n_runs"] >= MIN_RUNS_FOR_FIG].copy()
genes_with_key = set(agg_fig[agg_fig["contrast"] == "geno_trt"]["target"])
agg_fig = agg_fig[agg_fig["target"].isin(genes_with_key)]
agg_fig = agg_fig[~agg_fig["target"].isin(GENES_EXCLUDE)].copy()

# Per-run dots: filter to figure-subset genes
ddcq_fig = merged_per_run[
    merged_per_run["target"].isin(genes_with_key) &
    ~merged_per_run["target"].isin(GENES_EXCLUDE)
].copy()

print(f"Figure genes: {sorted(genes_with_key - set(GENES_EXCLUDE))}")
print(f"Rows in figure subset: {len(agg_fig)}")

make_figure(agg_fig, ddcq_fig,
            SUPP_FIG / "qPCR_figure_panel.pdf",
            SUPP_FIG / "qPCR_figure_panel.png")

# 10. Final summary
print("\n" + "=" * 72)
print("FINAL FIGURE-SUBSET SUMMARY (genes with >= 2 runs in geno_trt)")
print("=" * 72)
disp = agg_fig[["target", "contrast", "log2fc_mean", "log2fc_sem",
                "n_runs", "p_value", "padj_BH"]].copy()
disp["log2fc_mean"] = disp["log2fc_mean"].round(2)
disp["log2fc_sem"] = disp["log2fc_sem"].round(2)
disp["p_value"] = disp["p_value"].apply(
    lambda x: f"{x:.2g}" if np.isfinite(x) else "")
disp["padj_BH"] = disp["padj_BH"].apply(
    lambda x: f"{x:.2g}" if np.isfinite(x) else "")
print(disp.to_string(index=False))

# qPCR vs RNA-seq Spearman across all shared points
from scipy import stats as _stats
cmp_fig = comp_new[comp_new["gene"].isin(genes_with_key) &
                   ~comp_new["gene"].isin(GENES_EXCLUDE) &
                   comp_new["qpcr_log2fc"].notna() &
                   comp_new["rnaseq_log2fc"].notna()].copy()
rho, p = _stats.spearmanr(cmp_fig["qpcr_log2fc"], cmp_fig["rnaseq_log2fc"])
print(f"\nPanel C Spearman (figure subset): rho={rho:.3f}  "
      f"n={len(cmp_fig)}  p={p:.4g}")
