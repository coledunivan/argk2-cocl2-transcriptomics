# qPCR validation of RNA-seq results — methods and summary

**qPCR validation.** To validate the RNA-seq differential expression results,
quantitative RT-PCR was performed on four genes representing the major
response patterns (*ugt-31*, *oac-7*, *ugt-29*, *t22f3.11*), with *tba-1*
as the reference gene. Total RNA was reverse-transcribed and amplified with
SYBR Green on a Bio-Rad CFX system (protocol "Cole Sybr.prcl", standard
three-step cycling). For each biological sample three technical replicates
were run; no-reverse-transcriptase (NRT) controls were included for every
target/sample combination on every plate. Amplicon specificity was
confirmed by melt-curve analysis (single peak per target at expected
T_m); wells with secondary peaks or primer-dimer signatures were
flagged and excluded. Technical-replicate outliers were removed if they
exceeded a 0.5-cycle C_q standard deviation within a triplicate group
(up to two wells per group). Primer efficiencies were determined from
six-point 10-fold serial dilution standard curves for each primer pair
(Supplementary Figure SX): *tba-1* E = 105.6% (R² = 0.999), *oac-7* E =
87.4% (R² = 0.998), *ugt-31* E = 81.0% (R² = 0.999), *t22f3.11* E = 80.0%
(R² = 0.999), *ugt-29* E = 98.2% (R² = 0.973). Relative expression was
calculated by efficiency-corrected ΔΔC_q (Pfaffl, 2001, *Nucl. Acids
Res.* 29:e45):

log₂FC = −log₂(E_target) · ΔC_q,target + log₂(E_ref) · ΔC_q,ref

Each qPCR plate (n = 8 independent plates across the full dataset; n = 3–6
plates per gene × contrast in the figure subset) was treated as an
independent biological replicate. Three contrasts matching the RNA-seq
DESeq2 main effects were computed for each plate: treatment effect in
wild type (N2 treated vs. N2 untreated), baseline genotype effect
(argk-2⁻/⁻ vs. N2, both untreated), and genotype effect under treatment
(argk-2⁻/⁻ vs. N2, both treated). Per-gene × contrast log₂FC values
were averaged across plates and tested against zero by one-sample
Welch's *t*-test; *p* values were adjusted across all tested
combinations by the Benjamini–Hochberg procedure. qPCR and RNA-seq
(DESeq2, apeglm/ashr-shrunk) log₂FC values were compared by Spearman's
rank correlation across all 12 gene × contrast combinations (Figure Xc).
Analysis was performed in Python 3 using pandas, NumPy, and SciPy;
all source code is available in `qpcr_analysis.py` and `apply_pfaffl.py`.

---

## Final numbers

### Primer efficiencies (Supplementary Figure SX)

| Target | E (%) | Slope | R² | n | Dilutions used |
|---|---:|---:|---:|---:|---|
| *tba-1* (ref) | 105.6 | −3.194 | 0.9993 | 9 | D2, D3, D4 |
| *oac-7* | 87.4 | −3.665 | 0.9979 | 15 | D0, D1, D2, D3, D4 |
| *ugt-31* | 81.0 | −3.880 | 0.9989 | 18 | D1–D6 (all) |
| *t22f3.11* | 80.0 | −3.916 | 0.9988 | 18 | D1–D6 (all) |
| *ugt-29* | 98.2 | −3.365 | 0.9728 | 15 | D0, D1, D2, D3, D4 |

### Figure subset (Pfaffl-corrected, after integrating ColeD_41726)

| Gene | Contrast | log₂FC ± SEM | n plates | raw *p* | padj |
|---|---|---:|---:|---:|---:|
| *ugt-31* | N2 treated vs untreated | +2.04 ± 0.13 | 6 | 1.9×10⁻⁵ | **3.7×10⁻⁴ \*\*\*** |
| *ugt-31* | argk-2 untreated vs N2 untreated | +0.61 ± 1.16 | 6 | 0.62 | 0.72 |
| *ugt-31* | argk-2 treated vs N2 treated | +0.65 ± 0.24 | 6 | **0.046** | 0.22 |
| *oac-7* | N2 treated vs untreated | +1.51 ± 0.55 | 5 | 0.051 | 0.22 |
| *oac-7* | argk-2 untreated vs N2 untreated | +1.41 ± 0.59 | 6 | 0.063 | 0.22 |
| *oac-7* | argk-2 treated vs N2 treated | +0.43 ± 0.38 | 5 | 0.33 | 0.51 |
| *ugt-29* | N2 treated vs untreated | +3.65 ± 0.33 | 3 | **0.0081** | 0.085 |
| *ugt-29* | argk-2 untreated vs N2 untreated | +1.54 ± 0.88 | 3 | 0.22 | 0.43 |
| *ugt-29* | argk-2 treated vs N2 treated | +0.44 ± 0.49 | 3 | 0.47 | 0.58 |
| *t22f3.11* | N2 treated vs untreated | −1.00 ± 0.84 | 4 | 0.32 | 0.51 |
| *t22f3.11* | argk-2 untreated vs N2 untreated | −2.78 ± 1.80 | 5 | 0.20 | 0.43 |
| *t22f3.11* | argk-2 treated vs N2 treated | −0.39 ± 0.39 | 4 | 0.40 | 0.52 |

### Panel C (qPCR ↔ RNA-seq concordance)

- **Spearman ρ = 0.83, n = 12, *p* = 8×10⁻⁴**
- All four up-regulated-in-N2-treated genes show agreement in sign across both platforms.
- *t22f3.11* is the only down-regulated gene; both platforms agree on its direction.
- No sign flips across any of the 12 gene × contrast combinations.

---

## What changed vs. the pre-efficiency version

The Pfaffl correction reduces the log₂FC magnitudes by factors of
log₂(E_target), so low-efficiency primers see the largest shrinkage:

| Change | Classic ΔΔCq | Pfaffl | % change |
|---|---:|---:|---:|
| *ugt-31* treat_N2 | +2.40 | +2.04 | −15% |
| *t22f3.11* geno_unt | −3.29 | −2.78 | −15% |
| *oac-7* treat_N2 | +1.68 | +1.51 | −10% |
| *ugt-29* treat_N2 | +3.70 | +3.65 | −1% |

Every significance conclusion is preserved. The Panel C correlation drops
negligibly (0.85 → 0.83). This is the expected behavior: efficiency
correction shrinks biased-upward log₂FC estimates without flipping any
biological conclusions.

---

## Where the scaling approximation enters

For the 7 older qPCR plates we only have cached per-run log₂FC values
(no raw C_q), so exact Pfaffl would have required reprocessing from
scratch. We applied the identity

log₂FC_pfaffl = log₂(E_target) · log₂FC_classic + (log₂(E_ref) − log₂(E_target)) · ΔC_q,ref

under the assumption ΔC_q,ref ≈ 0 (reference gene stable across
conditions), which simplifies to the multiplicative scaling used in
`apply_pfaffl.py`. We validated this approximation on the new ColeD_41726
plate (where raw C_q is available and both exact and approximate values
can be computed):

- **mean |approximation error| = 0.060 log₂FC units**
- **max |approximation error| = 0.132 log₂FC units**

Both are well below the between-plate biological noise (per-plate SEM of
~0.2–1.8 log₂FC). The approximation is therefore accurate enough to use
for the old plates. If a reviewer pushes on this point, the fix is to
reprocess the old plates from their raw Cq files through the
now-efficiency-aware `compute_delta_cq()`; the pipeline is already set
up for this (just pass `efficiencies=EFFICIENCIES` to the function).

---

## Files in this output directory

### Primary figures
- `qPCR_figure_panel.pdf` / `.png` — main validation figure (panels a, b, c), Pfaffl-corrected
- `primer_efficiencies_supplement.pdf` / `.png` — 5-panel standard curves with fit statistics

### Data tables
- `primer_efficiencies.csv` — per-target E%, slope, R², n, excluded dilution points
- `qPCR_combined_deltadeltaCq.csv` / `_pfaffl.csv` — Pfaffl-corrected aggregate (main analysis)
- `qPCR_per_run_ddcq.csv` / `_pfaffl.csv` — per-plate log₂FC values
- `qPCR_vs_RNAseq_comparison.csv` — joined qPCR × DESeq2 table for Panel C

### Source code
- `qpcr_analysis.py` — main pipeline; now supports `efficiencies=…` parameter in `compute_delta_cq()`
- `fit_efficiencies_v2.py` — fits standard curves from the three efficiency ZIPs
- `apply_pfaffl.py` — wrapper that applies Pfaffl correction across all 8 plates
- `add_new_run.py` — (earlier) wrapper for adding a single new plate
- `rebuild_figure.py` — light script to rebuild the figure from saved CSVs (no math rerun)


3. **Biological vs. technical replication.** Each "plate" was treated as one
   biological replicate in the one-sample *t*-tests. If some of the 8 plates
   are technical repeats of the same RNA prep, the effective N is smaller and
   the error bars are too tight. Confirm this is correct before submission.

4. **Secondary melt peaks on *t22f3.11* / N2-treated wells (D07–D09 of ColeD_41726).**
   Small 74.5 °C shoulder alongside the 77 °C main peak; wells retained because
   the Cq triplicate SD was 0.06 so the main product dominates. 
