# Software versions used for the manuscript

These are the versions under which every figure in the manuscript was generated.
`install_R_packages.R` and `requirements.txt` install current releases, which
should be compatible, but this file is the authoritative record.

## R 4.5.1

| Package | Version | Role |
|---|---|---|
| DESeq2 | 1.42 | Differential expression, variance-stabilizing transformation |
| limma | — | Batch-effect removal for PCA visualization only |
| apeglm | — | LFC shrinkage (main coefficients) |
| ashr | — | LFC shrinkage (contrast combinations) |
| clusterProfiler | — | GO and KEGG enrichment |
| AnnotationDbi | — | Annotation database interface |
| org.Ce.eg.db | — | *C. elegans* gene annotation |
| celegans.db | — | *C. elegans* gene annotation |
| GO.db | — | Gene Ontology database |
| ggplot2 | — | Plotting |
| VennDiagram | — | Venn diagrams (Figure 1) |
| ggrepel | — | Label placement |
| pheatmap | — | Heatmaps |
| patchwork | — | Multi-panel figure assembly |
| RColorBrewer | — | Color palettes |
| scales | — | Axis formatting |
| tidyverse | — | dplyr, tidyr, readr, tibble, purrr, stringr |

## Python 3.10.14

| Package | Version | Role |
|---|---|---|
| matplotlib | 3.10.8 | Network rendering, figure assembly |
| numpy | 2.4.4 | Coordinate geometry, numerical operations |
| pandas | 3.0.2 | Node/edge table import and manipulation |
| scipy | — | Statistics for qPCR analysis |
| adjustText | — | Non-overlapping label placement |

---

**Recording exact versions.** Entries marked `—` were not pinned in the
manuscript's software section. To fill them in from the environment that
produced the figures:

```bash
# R
Rscript -e 'ip <- installed.packages()[, "Version"]; print(ip[c("limma","apeglm","ashr","clusterProfiler","AnnotationDbi","org.Ce.eg.db","celegans.db","GO.db","ggplot2","VennDiagram","ggrepel","pheatmap","patchwork","RColorBrewer","scales")])'

# Python
pip freeze | grep -iE "scipy|adjusttext"
```

A full `sessionInfo()` dump is printed at the end of
`environment/install_R_packages.R`; saving that output to
`environment/sessionInfo.txt` and committing it is good practice for a
reproducibility archive.
