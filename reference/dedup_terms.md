# Flag redundant terms for display

Related gene sets often reach significance together: in GO:BP or
Reactome, one biological signal can light up a dozen overlapping terms.
`dedup_terms()` picks one representative per cluster so a figure shows
the signal once. It changes what is drawn, not what was tested: no row
is dropped and no p-value is touched, because the multiple-testing
correction already covered every term.

## Usage

``` r
dedup_terms(
  enrichment,
  gene_sets,
  method = c("enrichmentmap", "collapse_pathways"),
  similarity = c("combined", "jaccard"),
  cutoff = NULL,
  term_threshold = 0.05,
  stats = NULL
)
```

## Arguments

- enrichment:

  An enrichment object from
  [`as_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_enrichment.md)
  or
  [`run_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/run_enrichment.md).

- gene_sets:

  A named list of character vectors, one per term, holding the full gene
  sets that were tested. Terms without a set are kept.

- method:

  `"enrichmentmap"` (gene-set overlap) or `"collapse_pathways"`
  (conditional enrichment).

- similarity:

  For `"enrichmentmap"`: `"combined"` or `"jaccard"`.

- cutoff:

  For `"enrichmentmap"`: similarity at or above which a term is
  redundant. `NULL` takes 0.375 for `"combined"` and 0.5 for
  `"jaccard"`.

- term_threshold:

  Only terms with `padj` below this are compared; the rest are left
  unflagged (`NA`). For `"collapse_pathways"` it is also the conditional
  p-value threshold.

- stats:

  For `"collapse_pathways"`: a named list with one ranking (named
  numeric vector of gene statistics) per contrast.

## Value

`enrichment`, with `dedup_status` (`"kept"`, `"redundant"` or `NA`),
`merged_into` and `similarity` columns in its results, and the settings
stored in `metadata$dedup`.

## How terms are compared

Within each contrast and database, significant terms are walked from the
smallest `padj` up. A term is `"redundant"` when its gene set is at
least `cutoff` similar to a representative already kept, and is recorded
as merged into the most similar one; otherwise it becomes a
representative. Similarity uses the full gene sets, not leading edges.

- `"combined"` (default, cutoff 0.375) is the EnrichmentMap coefficient:
  the mean of the Jaccard index and the overlap coefficient (Merico et
  al. 2010; Reimand et al. 2019). The overlap coefficient catches a
  small set nested inside a large one, which Jaccard misses.

- `"jaccard"` (cutoff 0.5) uses shared genes over all genes in either
  set.

## `method = "collapse_pathways"`

Runs
[`fgsea::collapsePathways()`](https://rdrr.io/pkg/fgsea/man/collapsePathways.html)
within each contrast and database: a term is redundant when it is no
longer enriched once the genes of a more significant term are
conditioned on. This is a statistical criterion, not an overlap rule, so
it needs the ranking each contrast was tested on and applies only to
ranked GSEA results (fgsea, clusterProfiler). It permutes, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible flags. The p-values it computes decide redundancy only; the
reported `padj` stay those of the original run.

## References

Merico D, Isserlin R, Stueker O, Emili A, Bader GD (2010). Enrichment
Map: a network-based method for gene-set enrichment visualization and
interpretation. PLoS ONE 5(11):e13984.
[doi:10.1371/journal.pone.0013984](https://doi.org/10.1371/journal.pone.0013984)

Reimand J, Isserlin R, Voisin V, et al. (2019). Pathway enrichment
analysis and visualization of omics data using g:Profiler, GSEA,
Cytoscape and EnrichmentMap. Nature Protocols 14:482-517.
[doi:10.1038/s41596-018-0103-9](https://doi.org/10.1038/s41596-018-0103-9)
