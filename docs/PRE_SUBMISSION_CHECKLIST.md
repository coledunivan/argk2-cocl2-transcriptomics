# Pre-submission checklist

Findings from checking the manuscript against the pipeline outputs in this repository.
Paragraph numbers in brackets refer to non-empty paragraphs in `Cole_Dunivan_G3.docx`,
counting from 0.

---

## First: what checks out

Worth knowing before the list of problems. Every headline statistic reproduces exactly
from the deposited data:

| Claim | Manuscript | Recomputed |
|---|---|---|
| G×E interaction genes | 90 | 90 ✓ |
| Sub-additive / supra-additive split | 68 / 22 | 68 / 22 ✓ |
| G×E genes in the TFLink network | 44 of 90 | 44 ✓ |
| G×E genes with ≥1 TF edge | 81 of 90 | 81 ✓ |
| Hub TFs rendered | 13 | 13 ✓ |
| qPCR ↔ RNA-seq concordance | ρ = 0.83, p = 0.0008, n = 12 | ρ = 0.832, p = 0.00079 ✓ |
| Primer efficiencies (5 targets) | 105.6 / 87.4 / 81.0 / 80.0 / 98.2 % | all ✓ |
| Primer R² | 0.9993 / 0.9979 / 0.9989 / 0.9988 / 0.9728 | all ✓ |

The analysis is sound and reproducible. Everything below is presentation, completeness,
and two genuine numerical discrepancies.

---

## Blockers — fix before your PI reads it

**1. Figures are not in the document.** Four placeholders are still sitting in Results:

- `[102]` `--Insert new figure 1 here`
- `[106]` `Insert new figure 2 here ---`
- `[110]` `Insert new figure 3 here --`
- `[115]` `Insert new figure 4 here ----`

Figures 5 and 6 have captions but no placeholder and no figure. Embed all six from
`outputs/figures/Figure1..6/` — use the `.tiff` at 600 dpi where available, otherwise
the `.pdf`.

**2. Missing required sections.** G3 requires these and the manuscript has none of them:

- Funding statement (grant numbers, or an explicit "no external funding" declaration)
- Author Contributions (CRediT format)
- Competing Interests statement
- Corresponding author footnote — paragraph `[2]` has `Cole Dunivan1,*` but the asterisk
  never resolves to an email address
- Keywords (5–8, after the Abstract)
- Acknowledgments — the CGC must be acknowledged; it is NIH-funded (P40 OD010440) and
  acknowledgment is a condition of strain distribution. You use three CGC strains.

**3. Data Availability still says "available upon request."** Replacement text with
GSE333535 already wired in is in `docs/DATA_AVAILABILITY.md`.

**4. First person singular in a two-author paper.** Three places:

- `[14]` "This allowed **me** to distinguish between the transcriptomic differences…"
- `[14]` "It also allowed **me** to parse non-additive gene-by-environment interactions…"
- `[105]` "we proceeded only with argk-2 for the rest of **my** study"

Change all to "us" / "our". `[105]` is the worst — it mixes both voices in one sentence.

**5. Two sentences begin with lowercase "we".**

- `[11]` "…and its cytosolic counterpart. **we** sought to determine…"
- `[14]` "…treatment or untreated conditions. **we** used a DESeq2 model…"

**6. "SYBER" should be "SYBR".** Appears twice — Methods `[91]` and the Supplementary
Figure 1 caption `[168]`.

**7. Doubled words.**

- `[168]` "qPCR-derived **derived** log₂(fold change)"
- `[170]` "**the the**"
- `[171]` "**expectation expectation**"

**8. Missing word.** `[105]` "To gain further into the molecular pathways" — presumably
"to gain further **insight** into".

---

## Numerical discrepancies

**9. Gene count is off by one.** `[93]` says "90 interaction genes from the **13,537**
that passed threshold filters." `outputs/DESeq2/DESeq2_effect_classification.csv` has
**13,536** rows. Change it, or say where the extra gene went.

**10. The "2,077 additive genes" in the Supplementary Figure 4 caption `[171]` is
mislabeled.** The effect classification is:

| Class | n |
|---|---|
| NS | 11,261 |
| Treatment | 2,013 |
| Genotype | 108 |
| G×E | 90 |
| Additive | 64 |

2,077 = Treatment (2,013) + Additive (64), which is exactly what
`12_Supp_additive_specificity.py` plots. But the caption calls all 2,077 "genes showing
significant additive" effects, and only **64** are classified Additive. Suggested
rewording: *"…for the 2,077 genes classified as Treatment-responsive or Additive…"*
As written, a reviewer who opens the CSV will think the number is wrong.

**11. The sample count does not add up, and the argk-1 exclusion is never stated.**

- Methods `[26]`/`[28]`: 32 libraries sequenced, 30 passed FastQC
- `data/Sample_Metadata_Table.csv`: **24** samples (3 genotypes × 2 treatments × 2
  experiments × 2 replicates)

32 samples is 4 genotypes × 8. The 8 argk-1 libraries were dropped, which the code
comments confirm ("argk-1 excluded — incomplete knockdown in raw counts") but the
manuscript never says. Paragraph `[17]` lists argk-1 among the strains grown, then says
"The RNA-seq analysis was performed on N2, RB2060, and RB2598" without explaining that
argk-1 *was* sequenced and then excluded on QC grounds.

This matters more than it looks: **GEO GSE333535 will contain the argk-1 libraries**, so
a reviewer comparing the accession to the manuscript will find eight libraries that the
paper does not account for. Add one sentence to Methods stating that argk-1 was
sequenced, showed incomplete transcript loss in raw counts, and was excluded before
modeling.

**12. Strain identifier conflict for argk-1.** Manuscript `[17]` says
`argk-1(ok2993) (MAH205)`; the header comment in `scripts/02_Figure1.R` says `MAH172`.
One is wrong. Worth resolving even though the strain is excluded, because it appears in
both the paper and the deposited code.

---

## Figures and supplementary material

**13. Supplementary Figures 3, 4 and 5 are never cited in the body text.** Each appears
exactly once — in its own caption. Journals require every supplementary item to be called
out in the main text. Natural homes:

- **Supp Fig 3** (90 G×E heatmaps) → the G×E paragraph `[133]`
- **Supp Fig 4** (sub-additive specificity) → `[133]`, where you argue the collapse is
  G×E-specific rather than global. This figure is the direct evidence for that claim and
  currently goes uncited.
- **Supp Fig 5** (dispersion estimates) → Methods `[57]`, which describes exactly this plot

**14. Supplementary Figure 3 has no script in the repository.** Supp 1 → `scripts/10`,
Supp 2 → `scripts/09`, Supp 4 → `scripts/12`, Supp 5 → `scripts/00`. Nothing generates
the 90 G×E heatmaps. Either the script was lost or the figure was made by hand — either
way it is the one panel a reader cannot reproduce from the deposit.

**15. `outputs/QC/` does not exist yet.** `00_QC_report.R` has never been run in this
tree, so Supplementary Figure 5 (dispersion) and the mapping-rate table are not in the
repository. Run it.

**16. Orphan figure file.** `outputs/figures/Figure1/Figure1_PanelG_Venn_Treatment.pdf`
is byte-identical to `Figure1_PanelD_Venn_Treatment.pdf` — a leftover from when the Venn
was panel G. The current script only writes panel D. Delete the panel G file so the
deposit does not ship a figure with no caption.

**17. `scripts/13_Supp_GxE_core_analysis.py` writes files named
`Figure6_panelB_GxE_indegree` and `Figure6_panelC_TF_coregulation`, but Figure 6's actual
panels B and C are the curated regulatory framework and the model schematic.** The names
collide with real panels that are different figures. This analysis (in-degree
distribution, 13 × 13 TF Jaccard heatmap) is not cited anywhere in the manuscript. Either
promote it to a numbered supplementary figure and cite it, or rename the outputs so
nobody mistakes them for Figure 6 panels.

**18. Stale TODO comment in `scripts/02_Figure1.R`.** The header block says *"DEGs for
Panel D: BH-adjusted p < 0.05 (no LFC cutoff)… Update figure caption to match."* The code
actually applies both `padj < 0.05` and `|log2FC| > 1` (line 608 onward), which matches
the caption. The caption is right and the comment is out of date — but a reviewer reading
the deposited code will hit that comment and think the caption overstates the filter.
Delete it.

---

## Also worth a look

- **Methods `[30]`–`[35]`, the HISAT2 command.** Line `[33]` uses
  `/pathtodirectory/cocl2_rnaseq/01.RawData/...` while `[34]` uses
  `/pathtodirectory/01.RawData/...` — inconsistent paths in what is presented as one
  command. Also it hardcodes sample A30; consider showing it generically.
- **Package versions.** Methods pins DESeq2, matplotlib, numpy and pandas but not limma,
  clusterProfiler, apeglm, ashr or scipy. `environment/versions.md` has a one-liner that
  pulls exact versions.
- **Figure 6 panel C is a hand-drawn schematic**, not script-generated. That is fine, but
  say so in the caption or legend so nobody looks for it in the code.
- **`[112]` is a Results subheading formatted as body text**, unlike the other subheadings
  (`[101]`, `[104]`, `[108]`, `[117]`). Cosmetic, but obvious in a Word file.

---

## Suggested order of work

1. Run `Rscript scripts/00_QC_report.R` → populates `outputs/QC/` (fixes #15)
2. Fix the text issues — #4 through #8, #9, #10 — one editing pass
3. Add the missing sections (#2) and swap in the new Data Availability text (#3)
4. Write the argk-1 exclusion sentence and resolve the strain ID (#11, #12)
5. Add the three supplementary citations (#13)
6. Embed all six figures (#1)
7. Clean the repo: delete the orphan panel G, fix the stale comment, decide on
   script 13's outputs (#16, #17, #18)
8. Decide about Supplementary Figure 3 (#14) — recover the script or note it as
   manually assembled
9. Push to GitHub and mint the Zenodo DOI (`docs/ZENODO.md`), then paste the DOI into
   the Data Availability statement

Items 1–6 are what your PI will notice. Items 7–9 are what a reviewer will notice.
