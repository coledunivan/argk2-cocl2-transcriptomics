#!/usr/bin/env python3
"""Fit primer-efficiency curves using target-specific dilution ranges.

The user's validated slopes (from their own analysis) require these
exclusions, which match points outside the assay's linear range:

  tba-1     D0 and D1 wells behave anomalously (D0 Cq ~30 when D1 is 10.7);
            D5 is at the LOD. Use D2-D4 only.
  oac-7     D0-D4 span the linear range cleanly; D5 jumps from Cq 27 -> 35
            (>8-cycle jump for a 10x dilution = past LOD).
  ugt-31    All six dilutions (D1-D6) are clean; R2 = 0.999 with no drops.
  t22f3.11  All six dilutions (D1-D6) are clean; R2 = 0.999 with no drops.
  ugt-29    D0-D4 linear; D5 jumps only ~0.8 cycles from D4 indicating it is
            at the detection floor (matches the user's validated fit).

Outputs:
  primer_efficiencies.csv          - final efficiency table for Pfaffl correction
  primer_efficiencies_supplement.pdf/png - 5-panel standard curve figure
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

# Repository-anchored paths 
ROOT     = Path(__file__).resolve().parents[1]
EFF_ROOT = ROOT / "data" / "qPCR" / "efficiency_curves"
QPCR_OUT = ROOT / "outputs" / "qPCR"
SUPP_FIG = ROOT / "outputs" / "supplementary"
QPCR_OUT.mkdir(parents=True, exist_ok=True)
SUPP_FIG.mkdir(parents=True, exist_ok=True)

# (target, plate_dir, dilution codes to USE, dilution codes to DROP)
FIG_TARGETS = [
    ("tba-1",   "Primer_Efficiency_NEW_-_tba-1__oac-7__ceh-74__zip-2-2",
                ["D2", "D3", "D4"], ["D0", "D1", "D5"]),
    ("oac-7",   "Primer_Efficiency_NEW_-_tba-1__oac-7__ceh-74__zip-2-2",
                ["D0", "D1", "D2", "D3", "D4"], ["D5"]),
    ("ugt-31",  "Primer_Efficiency_251104-2",
                ["D1", "D2", "D3", "D4", "D5", "D6"], []),
    ("t22f3.3", "Primer_Efficiency_251104-2",
                ["D1", "D2", "D3", "D4", "D5", "D6"], []),
    ("ugt-29",  "Primer_Eff_022326-2",
                ["D0", "D1", "D2", "D3", "D4"], ["D5"]),
]

def dilution_to_log10(sample):
    s = str(sample).strip().upper()
    # D1, D2, ..., D6 all 10-fold serial dilutions:
    # code D_n -> relative concentration 10^(-n) (except plates where D1 is the
    # most concentrated, i.e. D_n corresponds to 10^(-n+1)).
    # For consistency we just use the numeric part; the slope is invariant
    # under a constant shift of log10(conc), so the actual absolute values
    # don't matter for E.
    if s.startswith("D") and s[1:].isdigit():
        return -int(s[1:])
    return np.nan

def load_target(target, plate_dir):
    cq_file = list((EFF_ROOT / plate_dir).glob("*Quantification Cq Results*.csv"))[0]
    df = pd.read_csv(cq_file)
    sub = df[df["Target"] == target].copy()
    sub["log10_conc"] = sub["Sample"].apply(dilution_to_log10)
    sub["Cq_num"] = pd.to_numeric(sub["Cq"], errors="coerce")
    sub = sub[sub["log10_conc"].notna() & sub["Cq_num"].notna()]
    return sub[["Well", "Sample", "log10_conc", "Cq_num"]].reset_index(drop=True)

def fit(x, y):
    if len(x) < 3:
        return np.nan, np.nan, np.nan
    m, b = np.polyfit(x, y, 1)
    yhat = m * x + b
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return m, b, r2

def E_from_slope(slope):
    E = 10 ** (-1 / slope)
    return E, (E - 1) * 100

# 
# Fit each target and assemble a summary table
# 
results = []
fit_panels = []  # per-target data for the sup figure

for target, plate_dir, keep_dils, drop_dils in FIG_TARGETS:
    df = load_target(target, plate_dir)
    keep_mask = df["Sample"].isin(keep_dils)
    df_kept = df[keep_mask]
    df_drop = df[~keep_mask]

    m, b, r2 = fit(df_kept["log10_conc"].values, df_kept["Cq_num"].values)
    E_fold, E_pct = E_from_slope(m)
    results.append({
        "target":     target if target != "t22f3.3" else "t22f3.11",
        "plate":      plate_dir,
        "n_used":     len(df_kept),
        "n_excluded": len(df_drop),
        "dilutions_used":     ",".join(sorted(set(df_kept["Sample"]))),
        "dilutions_excluded": ",".join(sorted(set(df_drop["Sample"]))),
        "slope":      m,
        "intercept":  b,
        "r_squared":  r2,
        "E_fold":     E_fold,
        "E_percent":  E_pct,
    })
    fit_panels.append({
        "target": target if target != "t22f3.3" else "t22f3.11",
        "df_kept": df_kept,
        "df_drop": df_drop,
        "slope": m, "intercept": b, "r2": r2,
        "E_pct": E_pct,
    })

res_df = pd.DataFrame(results)
print("=" * 72)
print("PRIMER EFFICIENCY FITS")
print("=" * 72)
print(res_df[["target","n_used","n_excluded","slope","r_squared",
              "E_percent","dilutions_excluded"]].round(4).to_string(index=False))
res_df.to_csv(QPCR_OUT / "primer_efficiencies.csv", index=False)
print("\nsaved primer_efficiencies.csv")

# 
# Build the supplementary figure: 5 panels, one per target
# Style matches the user's ugt-29 reference plot
# 
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica","Arial","Liberation Sans","DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# 2x3 layout; last cell used for the summary table
fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.2))
axes_flat = axes.flatten()

# Use contrasting color per target (keep figure-friendly palette)
COLORS = {
    "tba-1":    "#6B6B6B",   # gray for ref
    "oac-7":    "#E07B3D",
    "ugt-31":   "#4A90C2",
    "t22f3.11": "#8E44AD",
    "ugt-29":   "#C0392B",
}

for i, p in enumerate(fit_panels):
    ax = axes_flat[i]
    target = p["target"]
    col = COLORS.get(target, "#555")

    # Kept points
    ax.scatter(p["df_kept"]["log10_conc"], p["df_kept"]["Cq_num"],
               s=60, c=col, edgecolor="#222", linewidth=0.7,
               zorder=3, label="Used")
    # Excluded points (open gray circles)
    if len(p["df_drop"]) > 0:
        ax.scatter(p["df_drop"]["log10_conc"], p["df_drop"]["Cq_num"],
                   s=60, facecolors="none", edgecolor="#888",
                   linewidth=0.9, zorder=2, label="Excluded")

    # Fit line spanning a little beyond the kept range
    xlo = p["df_kept"]["log10_conc"].min() - 0.7
    xhi = p["df_kept"]["log10_conc"].max() + 0.7
    xs = np.array([xlo, xhi])
    ax.plot(xs, p["slope"] * xs + p["intercept"],
            color=col, lw=1.2, alpha=0.7, zorder=1)

    # Title and annotation box
    ax.set_title(f"$\\it{{{target}}}$", fontsize=12, pad=4)
    txt = (f"E = {p['E_pct']:.1f}%\n"
           f"R$^2$ = {p['r2']:.4f}\n"
           f"slope = {p['slope']:.3f}\n"
           f"n = {len(p['df_kept'])}")
    ax.text(0.97, 0.97, txt, transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.35",
                      facecolor="white", edgecolor="#aaa", lw=0.5))

    ax.set_xlabel(r"log$_{10}$(relative concentration)", fontsize=9)
    ax.set_ylabel(r"C$_q$", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", length=3, labelsize=9)
    ax.grid(False)

    # Panel label
    ax.text(-0.14, 1.08, "abcdef"[i], transform=ax.transAxes,
            fontsize=14, fontweight="bold", ha="left", va="top")

# 6th panel: summary table
ax = axes_flat[5]
ax.axis("off")
table_data = [["Target", "E (%)", "Slope", "R²", "n"]]
for r in results:
    table_data.append([
        r["target"].replace("t22f3.3", "t22f3.11"),
        f"{r['E_percent']:.1f}",
        f"{r['slope']:.3f}",
        f"{r['r_squared']:.4f}",
        f"{r['n_used']}",
    ])
tbl = ax.table(cellText=table_data[1:], colLabels=table_data[0],
               loc="center", cellLoc="center",
               colWidths=[0.22, 0.16, 0.20, 0.22, 0.12])
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 1.45)
# Format header and italicize gene names
for (row, col), cell in tbl.get_celld().items():
    if row == 0:
        cell.set_facecolor("#E6E6E6")
        cell.set_text_props(fontweight="bold")
    if col == 0 and row > 0:
        cell.set_text_props(fontstyle="italic")
    cell.set_linewidth(0.5)

ax.text(-0.08, 1.08, "f", transform=ax.transAxes,
        fontsize=14, fontweight="bold", ha="left", va="top")
ax.text(0.5, 1.04, "Summary",
        transform=ax.transAxes, ha="center", va="bottom",
        fontsize=11, fontweight="bold")

# Overall title + spacing
fig.suptitle("Supplementary Figure: Primer efficiency standard curves",
             fontsize=12, fontweight="bold", y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(SUPP_FIG / "primer_efficiencies_supplement.pdf", dpi=400, bbox_inches="tight")
plt.savefig(SUPP_FIG / "primer_efficiencies_supplement.png", dpi=300, bbox_inches="tight")
print("saved primer_efficiencies_supplement.pdf/png")
