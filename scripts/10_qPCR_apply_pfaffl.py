#!/usr/bin/env python3
"""

Strategy

Primer efficiency correction modifies the per-well normalized expression:

  eff_dcq = Cq_target * log2(E_target) - mean(Cq_ref) * log2(E_ref)

which gives log2FC = -(mean(eff_dcq)_num - mean(eff_dcq)_den)
                  = -log2(E_target) * DeltaCq_target + log2(E_ref) * DeltaCq_ref
                  = Pfaffl (2001, NAR 29:e45).



  log2fc_pfaffl = log2(E_target) * log2fc_classic
                + (log2(E_ref) - log2(E_target)) * ΔCq_ref


"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "qpcr"))
from qpcr_analysis import (
    load_all_cq, flag_outliers, compute_delta_cq, compute_ddcq_contrasts,
    aggregate_across_runs, bh_adjust, make_figure,
    TARGET_RENAMES, REF_GENE, MIN_RUNS_FOR_FIG, GENES_EXCLUDE,
)

NEW_RUN = "ColeD_41726"
from qpcr_analysis import QPCR_OUT, QPCR_ROOT, SUPP_FIG  # repo-anchored paths

PER_RUN_SOURCE_CSV  = QPCR_OUT / "qPCR_per_run_ddcq_classic.csv"   # input  (never written)
PER_RUN_CLASSIC_CSV = QPCR_OUT / "qPCR_per_run_ddcq.csv"           # output (overwritten)
COMBINED_CSV        = QPCR_OUT / "qPCR_combined_deltadeltaCq.csv"
COMPARISON_CSV      = QPCR_OUT / "qPCR_vs_RNAseq_comparison.csv"
EFF_CSV             = QPCR_OUT / "primer_efficiencies.csv"

# 
# 0. Load primer efficiencies
# 
eff_df = pd.read_csv(EFF_CSV)
EFFICIENCIES = dict(zip(eff_df["target"], eff_df["E_fold"]))
# Fill in a default of 2.0 (100% efficient) for targets with no curve data.
# Among figure targets, all five have curves. 

for k, v in EFFICIENCIES.items():
    print(f"  {k:10s}  E={v:.3f} fold  ({(v-1)*100:.1f}%)  log2(E)={np.log2(v):.3f}")

REF_LOG2E = np.log2(EFFICIENCIES[REF_GENE])
print(f"\nRef log2(E) = {REF_LOG2E:.3f}  (classic assumes 1.000)")

# 1. Reprocess NEW run with exact Pfaffl correction
print("\n" + "=" * 72)
print("STEP 1: exact Pfaffl on new run (raw Cq available)")
print("=" * 72)

cq_df = load_all_cq(QPCR_ROOT)
assert cq_df["run"].nunique() == 1 and cq_df["run"].iloc[0] == NEW_RUN, \
    f"Expected only new run in workdir, got {cq_df['run'].unique()}"
cq_df = flag_outliers(cq_df)
n_dropped = (~cq_df["keep"]).sum()
print(f"[QC] {n_dropped} technical-replicate outlier well(s) dropped")

# Classic (reproduces the previous analysis) AND Pfaffl (new)
dcq_classic = compute_delta_cq(cq_df, ref_gene=REF_GENE, efficiencies=None)
dcq_pfaffl  = compute_delta_cq(cq_df, ref_gene=REF_GENE, efficiencies=EFFICIENCIES)
ddcq_classic = compute_ddcq_contrasts(dcq_classic)
ddcq_pfaffl  = compute_ddcq_contrasts(dcq_pfaffl)
for df in (ddcq_classic, ddcq_pfaffl):
    df["target"] = df["target"].replace(TARGET_RENAMES)

# Keep only the summary columns
classic_new = ddcq_classic[["run","target","contrast","ddcq","log2fc",
                            "n_num","n_den"]].copy()
pfaffl_new  = ddcq_pfaffl[["run","target","contrast","ddcq","log2fc",
                           "n_num","n_den"]].copy()

print("\nNew run log2FCs - classic vs exact Pfaffl:")
merged = classic_new.merge(pfaffl_new,
                           on=["run","target","contrast","n_num","n_den"],
                           suffixes=("_classic","_pfaffl"))
merged["delta"] = merged["log2fc_pfaffl"] - merged["log2fc_classic"]
print(merged[["target","contrast","log2fc_classic","log2fc_pfaffl","delta"]]
      .round(3).to_string(index=False))

print("\n" + "=" * 72)
print("STEP 2: scaling approximation on old runs (no raw Cq available)")
print("=" * 72)

if not PER_RUN_SOURCE_CSV.exists():
    sys.exit(
        f"[abort] {PER_RUN_SOURCE_CSV} not found.\n"
        "        This file holds the uncorrected per-plate DeltaDeltaCq values and is\n"
        "        the input to efficiency correction. It ships with the repository."
    )
old_per_run = pd.read_csv(PER_RUN_SOURCE_CSV)

old_per_run = old_per_run[old_per_run["run"] != NEW_RUN].copy()
print(f"  Old per-run rows: {len(old_per_run)}  (runs: {old_per_run['run'].nunique()})")

# Apply the scaling approximation. When ΔCq_ref ≈ 0, log2fc_pfaffl ≈
# log2(E_target) * log2fc_classic.  Fall back to 2.0 for targets we have
# no efficiency for (those genes won't appear in the main figure anyway).
def scale_factor(target):
    return np.log2(EFFICIENCIES.get(target, 2.0))

old_per_run["log2fc_classic"] = old_per_run["log2fc"]
old_per_run["log2fc"] = [
    scale_factor(t) * v
    for t, v in zip(old_per_run["target"], old_per_run["log2fc_classic"])
]
old_per_run["ddcq"] = -old_per_run["log2fc"]

print("\n  Validating scaling approximation vs exact Pfaffl on NEW run:")
val_rows = []
for _, r in classic_new.iterrows():
    exact = pfaffl_new[(pfaffl_new["target"] == r["target"]) &
                       (pfaffl_new["contrast"] == r["contrast"])]
    if exact.empty:
        continue
    approx = scale_factor(r["target"]) * r["log2fc"]
    val_rows.append({
        "target":   r["target"],
        "contrast": r["contrast"],
        "classic":  r["log2fc"],
        "approx":   approx,
        "exact":    exact["log2fc"].iloc[0],
        "approx_error": approx - exact["log2fc"].iloc[0],
    })
val_df = pd.DataFrame(val_rows)
print(val_df.round(3).to_string(index=False))
print(f"\n  approximation error: mean |err| = "
      f"{val_df['approx_error'].abs().mean():.3f}, "
      f"max |err| = {val_df['approx_error'].abs().max():.3f} log2FC units")

print("\n" + "=" * 72)
print("STEP 3: merging old (approx Pfaffl) + new (exact Pfaffl)")
print("=" * 72)

# Prepare old for merge
old_out = old_per_run[["run","target","contrast","ddcq","log2fc",
                       "n_num","n_den"]].copy()

combined = pd.concat([old_out, pfaffl_new], ignore_index=True)
print(f"  combined per-run rows: {len(combined)}  "
      f"({combined['run'].nunique()} plates)")

# 4. Re-aggregate and write CSVs
agg = aggregate_across_runs(combined)
agg["padj_BH"] = bh_adjust(agg["p_value"].values)
agg = agg.sort_values(["contrast","target"]).reset_index(drop=True)

# Save Pfaffl versions under distinct filenames so users can compare
combined.to_csv(QPCR_OUT / "qPCR_per_run_ddcq_pfaffl.csv", index=False)
agg.to_csv(QPCR_OUT / "qPCR_combined_deltadeltaCq_pfaffl.csv", index=False)
print(f"  wrote qPCR_per_run_ddcq_pfaffl.csv        ({len(combined)} rows)")
print(f"  wrote qPCR_combined_deltadeltaCq_pfaffl.csv ({len(agg)} rows)")

# Also overwrite the canonical names so downstream tools (figure builder,
# R comparison) pick up the corrected values by default
combined.to_csv(PER_RUN_CLASSIC_CSV, index=False)
agg.to_csv(COMBINED_CSV, index=False)
print(f"  overwrote {PER_RUN_CLASSIC_CSV} and {COMBINED_CSV} with Pfaffl-corrected values")

# 5. Refresh qPCR_vs_RNAseq_comparison.csv with Pfaffl-corrected qPCR columns
comp = pd.read_csv(COMPARISON_CSV)
# Collapse old duplicate rows (the t22f3.3 -> t22f3.11 rename had left
# a stale duplicate in the prior table)
comp = comp.drop_duplicates(subset=["gene","contrast"]).copy()

qpcr_map = agg.set_index(["target","contrast"]).to_dict("index")
def _pick(g, c, k):
    entry = qpcr_map.get((g, c))
    return entry.get(k) if entry else np.nan

comp["qpcr_log2fc"]  = [_pick(g, c, "log2fc_mean") for g, c in zip(comp["gene"], comp["contrast"])]
comp["qpcr_sem"]     = [_pick(g, c, "log2fc_sem")  for g, c in zip(comp["gene"], comp["contrast"])]
comp["qpcr_n_runs"]  = [_pick(g, c, "n_runs")      for g, c in zip(comp["gene"], comp["contrast"])]
comp["qpcr_pvalue"]  = [_pick(g, c, "p_value")     for g, c in zip(comp["gene"], comp["contrast"])]
comp["qpcr_padj_BH"] = [_pick(g, c, "padj_BH")     for g, c in zip(comp["gene"], comp["contrast"])]

sign_q = np.sign(comp["qpcr_log2fc"])
sign_r = np.sign(comp["rnaseq_log2fc"])
both = comp["qpcr_log2fc"].notna() & comp["rnaseq_log2fc"].notna()
comp["direction_agrees"] = np.where(both, sign_q == sign_r, pd.NA)
comp.to_csv(COMPARISON_CSV, index=False)
print(f"  refreshed {COMPARISON_CSV} ({len(comp)} rows)")

# 6. Figure-subset summary
print("\n" + "=" * 72)
print("FIGURE SUBSET (>= 2 runs per G x E, Pfaffl-corrected)")
print("=" * 72)
agg_fig = agg[agg["n_runs"] >= MIN_RUNS_FOR_FIG].copy()
keyc = "geno_trt"
genes_with_key = set(agg_fig[agg_fig["contrast"] == keyc]["target"])
agg_fig = agg_fig[agg_fig["target"].isin(genes_with_key) &
                  ~agg_fig["target"].isin(GENES_EXCLUDE)]
ddcq_fig = combined[combined["target"].isin(genes_with_key) &
                    ~combined["target"].isin(GENES_EXCLUDE)]

disp = agg_fig[["target","contrast","log2fc_mean","log2fc_sem",
                "n_runs","p_value","padj_BH"]].copy()
disp["log2fc_mean"] = disp["log2fc_mean"].round(2)
disp["log2fc_sem"] = disp["log2fc_sem"].round(2)
disp["p_value"] = disp["p_value"].apply(lambda x: f"{x:.2g}" if np.isfinite(x) else "")
disp["padj_BH"] = disp["padj_BH"].apply(lambda x: f"{x:.2g}" if np.isfinite(x) else "")
print(disp.to_string(index=False))

# 7. Build main figure
print("\n" + "=" * 72)
print("REBUILDING MAIN FIGURE (Pfaffl-corrected)")
print("=" * 72)
make_figure(agg_fig, ddcq_fig, "qPCR_figure_panel.pdf", "qPCR_figure_panel.png")

# 8. Spearman rho for panel C
from scipy import stats as _stats
cmp_fig = comp[comp["gene"].isin(genes_with_key) &
               ~comp["gene"].isin(GENES_EXCLUDE) &
               comp["qpcr_log2fc"].notna() &
               comp["rnaseq_log2fc"].notna()]
rho, p = _stats.spearmanr(cmp_fig["qpcr_log2fc"], cmp_fig["rnaseq_log2fc"])
print(f"\nPanel C Spearman rho = {rho:.3f}, n = {len(cmp_fig)}, p = {p:.4g}")
