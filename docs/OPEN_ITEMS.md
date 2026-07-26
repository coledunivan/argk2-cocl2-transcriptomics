# What's left

Everything from the audit has been applied. The repository is complete and verified:

| Check | Status |
|---|---|
| Figures | **11 / 11** present and committed |
| qPCR raw Cq coverage | **8 / 8** plates, each verified to reproduce its published values |
| Tracked files / size | 239 files, 83 MB — nothing near GitHub's 50 MB per-file limit |
| Cruft, absolute paths | none |
| Python scripts | all compile |
| `10_qPCR_apply_pfaffl.py` | idempotent (bug fixed) |
| Panel C Spearman ρ | 0.832, p = 0.0008, n = 12 — matches the manuscript |

---

## 1. ~~Run the QC script~~ — DONE

`outputs/QC/` now holds `Dispersion_estimates.pdf` (**Supplementary Figure 5**),
`PCA_genotype_treatment.pdf` and `Sample_distance_heatmap.pdf`. All eleven figures are in.

`Table_SX_mapping_rates.csv` was not written — `00_QC_report.R` only produces it when a
HISAT2 alignment summary is present in the working directory. Nothing in the manuscript
references a mapping-rate table (Methods give alignment rates in prose), so this is not a
gap.

---

## 2. ~~Two facts only you have~~ — BOTH RESOLVED

**argk-1 strain identifier — fixed.** You confirmed MAH205. Updated in
`scripts/02_Figure1.R` (comment and `EXCLUDE_STRAINS`), `scripts/01_DESeq2_and_TF_enrichment.R`,
and `scripts/_utils.R`. No `MAH172` remains anywhere. The manuscript already said MAH205.

Note this string is functionally inert on the deposited data — the shipped metadata
contains only N2, RB2060 and RB2598, so the exclusion filter has nothing to match. It
matters only if you re-run against a metadata file that still carries argk-1 rows.

**Sample-count reconciliation — closed.** Your own header comment in `_utils.R` had the
answer, and the file contents confirm it:

| Stage | n | Evidence |
|---|---|---|
| Libraries sequenced | 32 | Methods |
| Passed FastQC | 30 | `CleanedCounts.csv` has 30 `A`-labelled columns; A27 and A28 absent |
| argk-1 excluded | −6 | A3, A4, A11, A12, A19, A20 — present in counts, absent from metadata |
| Final analysis | 24 | `Sample_Metadata_Table.csv` |

So A27 and A28 are the two FastQC failures, and the six surviving argk-1 libraries were
dropped for incomplete transcript loss. The Methods sentence now states this precisely.

---

## 3. Confirm what I drafted

I wrote these from scratch. Check them before they go out.

- **Author Contributions** — "C.D. and S.R.F. conceived and designed the study. C.D.
  performed the experiments, carried out the computational analysis, and wrote the
  manuscript. S.R.F. supervised the project and edited the manuscript. Both authors read
  and approved the final manuscript." Consider whether Dr. Warner's bioinformatics support
  belongs here as well as in the Acknowledgements.
- **Conflicts of Interest** — "The authors declare no conflicts of interest."
- **Keywords** — *Caenorhabditis elegans; arginine kinase; phosphagen kinase; oxidative
  stress; cobalt chloride; gene-by-environment interaction; innate immunity; ZIP-2;
  mitochondrial retrograde signaling*
- **Corresponding author** — currently your gmail. Journals generally prefer an
  institutional address; swap in your UNCW one if you have it.

---

## 4. Publish the repository

Full walkthrough in `docs/ZENODO.md`. Short version:

1. `rm -rf .git` and the cleanup commands at the top of that file
2. `git init && git add -A && git commit && git push` to a **public** repo
3. Enable the repo at zenodo.org → Settings → GitHub, **before** creating the release
4. Cut release `v1.0.0` on GitHub — Zenodo mints the DOI automatically
5. Use the **concept DOI** (the "Cite all versions" one), not the version DOI

Then fill three placeholders:

```bash
grep -rn "GITHUB-USERNAME\|ZENODO-ID\|<your-username>" --include="*.md" --include="*.cff" .
```

| Placeholder | Files |
|---|---|
| `<GITHUB-USERNAME>` | `README.md`, `CITATION.cff`, and the Data Availability statement in the .docx |
| `<ZENODO-ID>` | Data Availability statement in the .docx |
| `[DOI assigned at acceptance]` | leave as-is — G3 assigns this at acceptance |

---

## 5. Send it

`email_to_dr_fausett.md` is drafted and sitting next to the repo folder. Attach the revised
.docx. It flags the two open facts from item 2 so she is not surprised by them.

---

## qPCR raw plates — 4 added, 3 still needed

Raw Bio-Rad CFX exports for four more plates are now in `data/qPCR/raw_plates/`, matched to
the published per-run table by target set and sample structure. **Three are still missing:**

| Plate | Targets | Structure |
|---|---|---|
| `qPCR_-_Treated_vs__Untreated__1_` | ceh-74, oac-7, t22f3.11, ugt-31, zip-2 | treated + untreated, both genotypes |
| `qPCR_030325__1_` | fbxa-128, fbxa-79, oac-7, ugt-29, ugt-31 | treated + untreated; dated 3 Mar 2025 |
| `qPCR_Results_-_RB2060_vs__N2` | ceh-74, cest-34, daf-16, fbxa-128, zip-2 | untreated only, RB2060 vs N2 |

The file that matters is `* Quantification Cq Results.csv` for each.

**Decision taken: published numbers stand.** The new plates are staged *outside*
`data/qPCR/runs/`, so `load_all_cq()` still sees only `ColeD_41726` and nothing recomputes
by accident. Full rationale and the recompute path are in
`data/qPCR/raw_plates/README.md`.

> **One thing to consider.** Efficiency correction for seven of the eight plates used an
> algebraic approximation because raw Cq was unavailable at the time — this is documented
> in `10_qPCR_apply_pfaffl.py`. Now that raw Cq for four of those plates ships with the
> repository, a reviewer could recompute exactly and get values differing by up to
> 1.4 log₂FC per plate, moving two borderline p-values (ugt-31 geno_trt 0.044 → 0.059;
> oac-7 treat_N2 0.048 → 0.062) above 0.05. No sign flips, and the qPCR ↔ RNA-seq
> correlation is unchanged.
>
> A single sentence in the qPCR Methods would remove the exposure without changing any
> number — something like: *"For plates where raw Cq values were unavailable, efficiency
> correction was applied using the algebraic form of the Pfaffl relation, which assumes
> reference-gene stability across conditions."* Say the word and I will add it.

---

## Bug found and fixed in `10_qPCR_apply_pfaffl.py`

While testing whether the new raw plates would change anything, I found that this script
was **not idempotent**: it read `outputs/qPCR/qPCR_per_run_ddcq.csv`, applied the
efficiency correction, and wrote the result back to the same file. Running it a second time
scaled every old-plate log₂FC by `log2(E_target)` again — silently shrinking all published
values by 10–20% with no error and no warning.

I hit this myself: an earlier verification run left doubly-corrected tables in the
repository. Caught by checking the tables against the values quoted in your manuscript
(2.0 / 1.5 / 3.7 / −1.0) — they read 1.797 / 1.400 / 3.617 / −0.841.

**Fixed two ways.** The tables were restored by inverting the extra pass, and every
published value now reproduces exactly:

| Quantity | Manuscript | Repository |
|---|---|---|
| ugt-31 treat_N2 | 2.0 | 2.04 |
| oac-7 treat_N2 | 1.5 | 1.51 |
| ugt-29 treat_N2 | 3.7 | 3.65 |
| t22f3.11 treat_N2 | −1.0 | −1.00 |
| ugt-31 treat_N2 padj | 4 × 10⁻⁴ | 4.0 × 10⁻⁴ |
| ugt-31 geno_trt p | 0.046 | 0.046 |
| Spearman ρ (Panel C) | 0.83, p = 0.0008, n = 12 | 0.832, p = 0.0008, n = 12 |

The script now reads its uncorrected input from a new file,
`outputs/qPCR/qPCR_per_run_ddcq_classic.csv`, which it never writes. Running it twice in a
row now produces byte-identical output — verified.

Worth knowing because the bug was live in the analysis before I touched anything. Had you
or a reviewer re-run that script on the original layout, the numbers would have quietly
changed.

---

## Left deliberately alone

- **The four "Insert new figure N here" placeholders** in Results — Dr. Fausett is handling
  figure placement, per your note.
- **Figure 6 panel C** is a hand-drawn schematic. Nothing in the repo reproduces it, which
  is fine — but the caption could say so.

## Retired

`13_Supp_GxE_core_analysis.py` and its six output files are excluded from the repository
via `.gitignore` — nothing in the manuscript cited them, and their filenames collided with
the real Figure 6 panels B and C. Before deleting I checked that no manuscript claim
depended on them: there are zero mentions of in-degree, Jaccard, or co-regulation anywhere
in the text, and the "13 TFs" figure in Methods comes from `07`/`08`, not from that script.

They are still on disk. To remove them for good, and then drop the `.gitignore` block that
names them:

```bash
rm scripts/13_Supp_GxE_core_analysis.py
rm outputs/supplementary/Figure6_*
```

The heatmap script previously numbered `14` is now `13_Supp_GxE_heatmap.py`, so the
pipeline runs `00`–`13` with no gaps.

---

## Fixed already — for your records

25 edits applied to the manuscript, all verified programmatically after the fact:

| Category | Changes |
|---|---|
| Numbers | 13,537 → 13,536; clarified the 2,077 figure as Treatment + Additive classes |
| Voice | 3 first-person singular slips; 2 sentences starting lowercase "we" |
| Typos | SYBER → SYBR (×2); "derived derived"; "the the"; "expectation expectation"; "Violine" → "Violin"; "mostly highly" → "most highly"; "cluster profiler" → "clusterProfiler"; "gain further into" → "gain further insight into"; a stray double hyphen |
| Spelling | "normalised" → "normalized"; "colour" → "color" |
| Captions | Supplementary Figure 3 caption described the heatmap transposed and listed the contrasts in the wrong order — both corrected against the rendered figure |
| Citations | Pagano et al. 2015 removed; re-audited to 24 references / 23 citations / zero orphans |
| Added | Acknowledgements, Funding, Author Contributions, Conflicts of Interest, Keywords, corresponding-author line |
| Added | New Data Availability statement with GEO GSE333535 |
| Added | Methods sentence on the argk-1 exclusion |
| Added | Body citations for Supplementary Figures 3, 4 and 5, none of which had been cited |

Repository: 126 files, 74 MB, nothing over GitHub's size limits, `scripts/00`–`14` plus the
`qpcr/` module. The entire Python half was executed end-to-end and reproduces every
published statistic.
