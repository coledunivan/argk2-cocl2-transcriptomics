# Data

Processed and reference data for the ARGK-2 / CoCl₂ RNA-seq study.
**Raw sequencing reads are not here** — they are deposited at NCBI GEO under
accession [GSE333535](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333535).

## Count matrices and metadata

| File | Description |
|---|---|
| `RNASEQ61125.csv` | Primary gene × sample count matrix. Rows are WormBase sequence names; columns are GEO library names (`<strain>_<treatment>_EXP<n>_R<n>`). Input to `01`, `02`, `04`. |
| `CleanedCounts.csv` | Annotated counts carrying `gene_symbol`, `ENTREZID`, `ENSEMBL`, `WORMBASE` and `GENENAME` alongside sample columns `A1`–`A32`. Retained for provenance and ID mapping. |
| `Sample_Metadata_Table.csv` | Sample sheet: `SampleLabel` (A1–A32), `genotype`, `treatment`, `experiment`, `replicate`, `GEO_library_name`. The `GEO_library_name` column is the join key between the two count matrices and the GEO submission. |

**Design.** 3 genotypes (N2, RB2060 = *argk-2(ok2723)*, RB2598 = *argk-4(ok3602)*)
× 2 treatments (untreated, 16 mM CoCl₂ for 3 h) × 2 independent experiments
× 2 biological replicates. Sample labels A27 and A28 are absent — those libraries
did not pass QC. `argk-1` (MAH205) was excluded from all analyses because the raw
counts showed incomplete knockdown.

## `reference/` — offline annotation resources

| File | Source | Used by |
|---|---|---|
| `RNAseq_to_TF_Targets.csv.gz` | [TFLink](https://tflink.net) TF–target pairs, joined to the count matrix. One row per TF–target edge, so genes repeat; scripts deduplicate on `Gene`. Gzipped (~6 MB; ~79 MB uncompressed). | `01`, `04`, `05`, `06`, `07` |
| `org_Ce_eg_GO_map.csv` | GO term → gene mapping exported from `org.Ce.eg.db` | `04`, `_utils.R` |
| `kegg_cel_pathway_to_gene.csv` | KEGG `cel` pathway → gene mapping | `03`, `_utils.R` |
| `kegg_cel_pathway_names.csv` | KEGG `cel` pathway ID → human-readable name | `03`, `_utils.R` |

These are committed so enrichment analyses reproduce exactly and do not depend on
live KEGG/GO web services, which change between releases.

> The `.gz` is read transparently by both R and pandas — no manual decompression
> needed. Scripts prefer an uncompressed `RNAseq_to_TF_Targets.csv` if one exists
> locally (it is gitignored), otherwise they use the `.gz`.

## `qPCR/` — Bio-Rad CFX exports

| Directory | Contents |
|---|---|
| `runs/ColeD_41726/` | One qPCR plate (all four genotype × treatment conditions): Quantification Cq Results, Melt Curve Summary, Melt Curve Peak Results, Run Information. |
| `efficiency_curves/Primer_Eff_022326-2/` | Standard curves: ugt-29, fbxa-79, fbxa-128, fbxc-58 (D0–D5) |
| `efficiency_curves/Primer_Efficiency_251104-2/` | Standard curves: ugt-31, t22f3.11, cest-34, fbxa-79 (D1–D6) |
| `efficiency_curves/Primer_Efficiency_NEW_-_tba-1__oac-7__ceh-74__zip-2-2/` | Standard curves: tba-1, oac-7, ceh-74, zip-2 (D0–D5) |

Raw Cq files for the seven earlier plates were not retained. Their per-plate ΔΔCq
values are preserved in `outputs/qPCR/qPCR_per_run_ddcq.csv`, and
`10_qPCR_apply_pfaffl.py` recomputes the aggregate across all eight plates from
that table plus the new run — see `docs/qPCR_integration_report.md` for the
consequences of this for exact reproducibility.

## License

Data in this directory are CC BY 4.0 (see `../LICENSE-DATA`). Redistributed
third-party resources retain their original terms, listed in that file.
