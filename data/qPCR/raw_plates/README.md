# Raw qPCR plates — complete set of all 8

Bio-Rad CFX exports for **all eight plates** behind Supplementary Figure 1. Directory names
match the run names in `../../../outputs/qPCR/qPCR_per_run_ddcq.csv` exactly.

`ColeD_41726` is also present in `../runs/`, which is what `load_all_cq()` scans and what
the published pipeline consumes. The other seven stay here rather than in `../runs/`
**on purpose** — moving them changes what the pipeline computes. See "Why this matters".


Two naming quirks worth knowing:

- **`qPCR_030325__1_` ran on 3 March 2026, not 2025.** The published run name reads `030325`
  but the CFX file is `ColeD_030326` and the instrument recorded 03/03/2026. The plate name
  carries a year typo; the run date above is authoritative.
- **`qPCR_-1_31_2026...` started on 30 January** and finished after midnight, which is where
  the `1_31_2026` in the published name comes from.
- **`qPCR_-_Treated_vs__Untreated__1_`** comes from a `.pcrd` named `Cole qPCR 250511`, which
  is the plate-setup template name, not the run date. Three other plates reference the same
  `Cole qPCR 250511.pltd` setup file.

Each directory holds the full CFX export set — Quantification Cq Results, Melt Curve
Summary / Peak Results / Derivative, Standard Curve Results, Run Information, and the
plate-view and end-point files. Run dates come from the `Run Information.csv` export in
each directory, so they are recorded in the deposit and need not be tracked separately.

Note the `ColeD_13026` plate was **started** on 30 January and **ended** after midnight on
31 January — which is where the published run name `qPCR_-1_31_2026...` comes from. That
resolves the only genuine ambiguity in the mapping, since this plate and
`ColeD_22026_qPCR` carry identical target sets.

### Assignment verified numerically

Directory names were assigned by recomputing classic ΔΔCq from each raw plate and comparing
against **every** published plate in `outputs/qPCR/qPCR_per_run_ddcq_classic.csv`, not by
filename. Result — an 8 × 8 cross-match with a unique 1:1 mapping:

| Raw plate | Best-matching published plate | max abs. difference | rows |
|---|---|---|---|
| `ColeD_41726` | `ColeD_41726` | 4.4 × 10⁻¹⁶ | 12 |
| `Fw__qPCR_-_022026` | `Fw__qPCR_-_022026` | 1.8 × 10⁻¹⁵ | 12 |
| `qPCR - Treated vs. Untreated` | `qPCR_-_Treated_vs__Untreated__1_` | 1.8 × 10⁻¹⁵ | 15 |
| `qPCR 030325` | `qPCR_030325__1_` | 1.8 × 10⁻¹⁵ | 15 |
| `qPCR Results - RB2060 vs. N2` | `qPCR_Results_-_RB2060_vs__N2` | 0.0 | 5 |
| `ColeD_13026` | `qPCR_-1_31_2026_-_tba-1...` | 4.4 × 10⁻¹⁶ | 12 |
| `Cole 251102` | `qPCR_-_tba-1__oac-7__ceh-74...` | 5.6 × 10⁻¹⁷ | 5 |
| `cole qPCR 040626` | `qPCR_4_6_26` | 1.3 × 10⁻¹⁵ | 9 |

Every plate reproduces its published counterpart to floating-point precision, and every
off-diagonal comparison is large (0.4–6.5 log₂FC). All 85 published per-plate values are
now backed by raw Cq.

## Coverage

All 8 of 8 published plates now have raw Cq data. Nothing is missing.

## Why this matters — read before moving anything into `runs/`

`10_qPCR_apply_pfaffl.py` documents that raw Cq was unavailable for seven of the eight
plates, so efficiency correction for those plates used an algebraic **approximation**:

```
log2fc_pfaffl ≈ log2(E_target) · log2fc_classic + (log2(E_ref) − log2(E_target)) · ΔCq_ref
```

which reduces to a simple multiplication only when the reference gene is perfectly stable
(ΔCq_ref ≈ 0). With raw Cq now available for four more plates, the exact Pfaffl values are
computable — and they differ from the published approximations:

| | |
|---|---|
| Rows compared | 50 |
| Max abs. difference | **1.39 log₂FC** |
| Mean abs. difference | 0.17 log₂FC |
| Rows differing > 0.05 | 31 of 50 |
| **Sign flips** | **0** |
| `ColeD_41726` (always had raw Cq) | 0.0000 — exact match, as expected |

**This gap is not caused by the three missing plates.** The classic (uncorrected) ΔΔCq
values computed from raw Cq match the published ones exactly — 0.000 on every diagonal
above. The divergence appears only at the efficiency-correction step, and it is a property
of the approximation itself.

The approximation assumes the reference gene is stable between the two conditions being
contrasted (ΔCq_ref ≈ 0). On these plates tba-1 moves by up to **2.53 cycles**:

| Plate | tba-1 mean Cq, N2 untreated → treated |
|---|---|
| `qPCR_4_6_26` | 22.51 → 25.03 (2.52 cycles) |
| `qPCR_-1_31_2026...` | 21.20 → 23.63 (2.43 cycles) |
| `Fw__qPCR_-_022026` | 26.77 → 26.08 (0.69 cycles) |
| `ColeD_41726` | 24.16 → 23.82 (0.34 cycles) |

Approximation error tracks that drift (correlation between |ΔCq_ref| and |error| = 0.67).
Obtaining the last three plates would let exact correction be applied to *those* plates
too; it would not shrink the error on these four.

Most of this drift is plausibly cDNA input differences rather than tba-1 being
treatment-responsive — ΔCq normalization within each sample is exactly what corrects for
loading — but it is the kind of thing a reviewer may ask about, and it is why the
approximation and the exact calculation part company.

Direction is preserved everywhere, so no biological conclusion reverses. But at the
aggregate level two significance calls in Supplementary Figure 1 sit right on the boundary
and move across it:

| Gene × contrast | Published | Recomputed | p (published) | p (recomputed) |
|---|---|---|---|---|
| ugt-31 geno_trt | 0.556 | 0.570 | 0.044 | **0.059** |
| oac-7 treat_N2 | 1.400 | 1.542 | 0.048 | **0.062** |

Both currently marked significant; both would fall just above 0.05. The manuscript cites
"ugt-31 geno_trt: p = 0.046" specifically.

Everything else is stable. `ugt-31 treat_N2` stays at p ≈ 4 × 10⁻⁴, `ugt-29 treat_N2` stays
significant, and the qPCR ↔ RNA-seq rank correlation is essentially unchanged (it moved
slightly *up* in a like-for-like comparison).

## The decision

**Option A — leave it.** Published numbers stand. The raw data still ships, so the deposit
is richer than before, and this README documents the approximation honestly. Nothing in the
manuscript changes.

**Option B — recompute once all eight plates are in hand.** Cleanest scientifically: every
plate gets exact Pfaffl correction and the whole of Supplementary Figure 1 becomes
reproducible from raw Cq. Costs two significance markers and requires updating the caption,
the p-values, and the qPCR paragraph in Results.

**Option C — recompute now with 5 exact + 3 approximated.** Not recommended. A mixed
aggregate is harder to describe in Methods than either clean option.

To do B once the last three plates arrive:

```bash
mv data/qPCR/raw_plates/* data/qPCR/runs/
python3 scripts/09_qPCR_fit_efficiencies.py
python3 scripts/10_qPCR_apply_pfaffl.py
```

and then update `10_qPCR_apply_pfaffl.py`, which currently branches on the assumption that
only one plate has raw Cq.
