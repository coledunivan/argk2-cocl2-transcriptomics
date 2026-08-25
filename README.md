# ARGK-2 sets a transcriptional stress response baseline in *Caenorhabditis elegans*

Analysis code and processed data for Dunivan & Fausett, submitted to *G3: Genes|Genomes|Genetics*.

Full-factorial bulk RNA-seq of *C. elegans* N2 (wild type), *argk-2(ok2723)* (RB2060) and
*argk-4(ok3602)* (RB2598) under acute CoCl₂ oxidative stress, with orthogonal qPCR validation.
Every figure in the manuscript is reproducible from this repository.

Resource | Location |
Raw sequencing reads | NCBI GEO [GSE333535](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333535) |
Code + processed data | This repository |
Archived snapshot (citable DOI) | Zenodo — 10.5281/zenodo.21613293 |

---

## Quick start

```bash
git clone https://github.com/coledunivan/argk2-cocl2-transcriptomics.git
cd argk2-cocl2-transcriptomics

# R dependencies (CRAN + Bioconductor)
Rscript environment/install_R_packages.R

# Python dependencies
pip install -r environment/requirements.txt

# Run the pipeline in order, always from the repository root
Rscript scripts/01_DESeq2_and_TF_enrichment.R
Rscript scripts/02_Figure1.R
# ... etc, see the table below
```

**All scripts must be launched from the repository root.** Each script resolves its own
project root, so `Rscript scripts/02_Figure1.R` and `cd scripts && Rscript 02_Figure1.R`
both work, but relative paths are always interpreted against the repo root.


## Pipeline
Scripts are numbered in dependency order. `01` must run first — it produces the DESeq2
contrast tables that almost everything downstream consumes.

| # | Script | Produces | Manuscript |
|---|---|---|---|
| 00 | `00_QC_report.R` | `outputs/QC/` — PCA, sample-distance heatmap, dispersion, mapping rates | QC (not shown) |
| 01 | `01_DESeq2_and_TF_enrichment.R` | `outputs/DESeq2/` — all five contrast tables, effect classification, TF target enrichment, `dds` object | Methods |
| 02 | `02_Figure1.R` | `outputs/figures/Figure1/` + per-contrast DEG tables | **Figure 1** |
| 03 | `03_Figure2_KEGG.R` | `outputs/figures/Figure2/` | **Figure 2** |
| 04 | `04_Figure3_GO.R` | `outputs/figures/Figure3/` | **Figure 3** |
| 05 | `05_Figure4_TF.R` | `outputs/figures/Figure4/` | **Figure 4** |
| 06 | `06_Figure5_GxE.py` | `outputs/figures/Figure5/` | **Figure 5** |
| 07 | `07_GRN_inputs.R` | `outputs/GRN_tables/` — nodes, edges, TFLink coverage report | Figure 6 input |
| 08 | `08_Figure6_GRN.py` | `outputs/figures/Figure6/` | **Figure 6** |
| 09 | `09_qPCR_fit_efficiencies.py` | `outputs/qPCR/primer_efficiencies.csv`, standard-curve panel | **Supp. Figure 2** |
| 10 | `10_qPCR_apply_pfaffl.py` | `outputs/qPCR/` — efficiency-corrected ΔΔCq tables, qPCR panel | **Supp. Figure 1** |
| 11 | `11_qPCR_rebuild_figure.py` | Re-renders the qPCR panel from saved CSVs (no recomputation) | — |
| 12 | `12_Supp_additive_specificity.py` | `outputs/supplementary/Figure_S_additive_vs_gxe.*` | Supp. figure |
| 13 | `13_Supp_GxE_heatmap.py` | `outputs/supplementary/Figure_S3_GxE_heatmap.*` | **Supp. Figure 3** |

Shared helpers live in `scripts/_utils.R` (sourced automatically) and `scripts/qpcr/`
Dependency order: `01` → everything. `07` → `08`. `09` → `10` → `11`.

Supplementary figure map: Supp 1 → `10`, Supp 2 → `09`, Supp 3 → `13`, Supp 4 → `12`,
Supp 5 → `00` (`outputs/QC/Dispersion_estimates.pdf`).

## Repository layout

```
argk2-cocl2-transcriptomics/
├── data/                              # inputs (see data/README.md)
│   ├── RNASEQ61125.csv                # raw gene × sample count matrix
│   ├── CleanedCounts.csv              # annotated counts (gene symbol, ENTREZ, ENSEMBL, WormBase)
│   ├── Sample_Metadata_Table.csv      # genotype × treatment × experiment × replicate
│   ├── reference/                     # offline annotation resources
│   │   ├── RNAseq_to_TF_Targets.csv.gz    # TFLink TF–target pairs joined to counts (gzipped, ~6 MB)
│   │   ├── org_Ce_eg_GO_map.csv           # GO term → gene map
│   │   ├── kegg_cel_pathway_to_gene.csv   # KEGG cel pathway → gene
│   │   └── kegg_cel_pathway_names.csv     # KEGG cel pathway IDs → names
│   └── qPCR/
│       ├── runs/ColeD_41726/          # Bio-Rad CFX exports (Cq, melt curve, run info)
│       └── efficiency_curves/         # three standard-curve plates
├── scripts/                           # numbered pipeline, see table above
│   ├── _utils.R                       # shared R helpers
│   └── qpcr/                          # qPCR library module + add_new_run.py utility
├── outputs/                           # everything the pipeline generates
│   ├── DESeq2/                        # contrast tables, effect classification, dds object
│   ├── GRN_tables/                    # network nodes / edges / coverage
│   ├── figures/Figure1..6/            # main figures + their source-data CSVs
│   ├── supplementary/                 # supplementary figures
│   └── qPCR/                          # efficiency and ΔΔCq tables
├── environment/                       # dependency manifests
├── docs/                              # data availability statement, Zenodo guide, qPCR report
├── LICENSE                            # MIT (code)
└── LICENSE-DATA                       # CC BY 4.0 (data and figures)
```

`outputs/` is committed so reviewers can inspect every figure and table without
re-running the pipeline. Re-running the scripts overwrites it in place.

## The gzipped TFLink file

`data/reference/RNAseq_to_TF_Targets.csv` is ~79 MB uncompressed, past the point where
GitHub is comfortable. The repository ships `RNAseq_to_TF_Targets.csv.gz` (~6 MB) instead.

No action is needed: R (`read.csv`, `readr::read_csv`) and Python (`pandas.read_csv`)
both read `.gz` transparently, and every script prefers an uncompressed copy if you have
one locally, falling back to the `.gz`. To decompress anyway:

```bash
gunzip -k data/reference/RNAseq_to_TF_Targets.csv.gz
```

---

## Requirements

**R ≥ 4.3** with Bioconductor 3.18+:

`DESeq2`, `apeglm`, `ashr`, `limma`, `clusterProfiler`, `org.Ce.eg.db`, `celegans.db`,
`GO.db`, `AnnotationDbi`, `tidyverse` (`dplyr`, `ggplot2`, `readr`, `stringr`, `tibble`,
`tidyr`, `purrr`), `patchwork`, `ggrepel`, `pheatmap`, `VennDiagram`, `RColorBrewer`, `scales`

**Python ≥ 3.10:** `pandas`, `numpy`, `scipy`, `matplotlib`, `adjustText`

`Rscript environment/install_R_packages.R` and `pip install -r environment/requirements.txt`
install both sets. Exact versions used for the manuscript are recorded in
`environment/versions.md`.

Script `07_GRN_inputs.R` downloads the WormBase gene-ID table (~2 MB) to map sequence
names to common names; it needs network access on first run.

## Citation

If you use this code or data, please cite the manuscript and the Zenodo archive.
See `CITATION.cff`, or use GitHub's "Cite this repository" button.

## License

Code (`scripts/`, `environment/`) is MIT — see `LICENSE`.
Data and figures (`data/`, `outputs/`) are CC BY 4.0 — see `LICENSE-DATA`.

## Contact

Sarah Fausett — Department of Biology and Marine Biology, University of North Carolina Wilmington
fausetts@uncw.edu
