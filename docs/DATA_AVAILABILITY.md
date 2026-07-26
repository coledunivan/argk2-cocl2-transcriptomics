# Data Availability statement

## What the manuscript currently says

> All R and Python scripts used for analysis are available upon request. Raw and
> processed count matrices, sample metadata, and differential expression results
> will be provided on GitHub.

Two problems for a G3 submission. "Available upon request" no longer satisfies most
journal data policies — G3 expects code and data to be deposited and accessible at
review. And "will be provided" is a future promise rather than a citable location;
reviewers need a working link now.

## Replacement text

Paste this over the existing Data Availability paragraph. The GitHub URL is already
filled in; substitute the Zenodo DOI once the release is archived.

---

**Data availability**

Raw sequencing reads and processed count matrices have been deposited in the NCBI
Gene Expression Omnibus (GEO) under accession GSE333535. All analysis code, processed
data, and figure source data are publicly available at
https://github.com/coledunivan/argk2-cocl2-transcriptomics and archived at Zenodo
(doi:10.5281/zenodo.XXXXXXX). The repository contains the complete numbered analysis
pipeline in R and Python, the annotated count matrices and sample metadata, the offline
Gene Ontology, KEGG and TFLink annotation resources used for enrichment analysis, all
DESeq2 contrast tables and effect classifications, the gene-regulatory-network node and
edge tables, the raw Bio-Rad CFX qPCR exports and primer-efficiency standard curves, and
every figure in this manuscript together with the source data underlying each panel.
Running the scripts in numbered order from the repository root reproduces all main and
supplementary figures. Strains are available from the Caenorhabditis Genetics Center.
Supplemental material available at figshare: [G3 will assign this DOI at acceptance].

---

## Notes on the wording

- **"Strains are available from the CGC"** — G3 requires a statement about strain and
  reagent availability, not just data. N2, RB2060 and RB2598 are CGC strains; the
  backcrossed FAU5 and FAU7 lines used for qPCR were generated in your lab, so if you
  want to be thorough, add: *"Backcrossed strains FAU5 and FAU7 are available from the
  corresponding author."*
- **The figshare line** is G3 boilerplate. G3 hosts supplemental material at figshare and
  assigns that DOI at acceptance — leave the bracketed placeholder in the submitted
  version.
- **Order matters.** GEO first (raw data), then GitHub (working code), then Zenodo
  (frozen citable archive). Reviewers look for the accession first.

## Before you submit — checklist

- [ ] GEO submission GSE333535 is released, or set to release on publication
- [ ] GitHub repository is **public**
- [ ] Zenodo DOI minted (see `docs/ZENODO.md`) and substituted for `10.5281/zenodo.XXXXXXX`
- [x] GitHub username filled in throughout (`coledunivan`)
- [ ] Methods section names the specific software versions — already present, see
      `environment/versions.md`
- [ ] The `argk-1`/MAH205 exclusion is stated in Methods (it is currently only a code comment)

## Also worth checking

The Methods section lists package versions for DESeq2, matplotlib, numpy and pandas but
not for limma, clusterProfiler, apeglm, ashr, or scipy. `environment/versions.md` has a
one-liner that pulls the exact versions from your R and Python environments — reviewers
of computational papers do sometimes ask.
