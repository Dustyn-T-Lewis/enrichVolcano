# Deduplicate enrichment results

Deduplicate enrichment results

## Usage

``` r
ev_collapse(
  enrich_result,
  method = c("jaccard", "collapse_fgsea", "both"),
  cutoff = 0.5,
  scope = c("within_db", "cross_db", "global"),
  collapse_pval_threshold = 0.05,
  keep_by = c("padj", "size", "NES")
)
```

## Arguments

- enrich_result:

  Tibble from
  [`ev_enrich()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_enrich.md).

- method:

  `"jaccard"`, `"collapse_fgsea"`, or `"both"`.

- cutoff:

  Jaccard similarity cutoff.

- scope:

  `"within_db"`, `"cross_db"`, or `"global"`.

- collapse_pval_threshold:

  For `collapse_fgsea`.

- keep_by:

  Tie-breaking column ("padj", "size", "NES").

## Value

Tibble with `dedup_kept` logical column (TRUE = retained).

## Details

- **Jaccard**: \\J(A,B) = \|A \cap B\| / \|A \cup B\|\\. Above cutoff,
  keep representative (by `keep_by`). Data-independent.

- **collapse_fgsea**:
  [`fgsea::collapsePathways`](https://rdrr.io/pkg/fgsea/man/collapsePathways.html)
  — tests leading-edge dependency. Data-dependent; preferred for fgsea
  output. Requires the `ev_pathways` attribute attached by
  [`ev_enrich()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_enrich.md);
  if absent, warns and keeps all pathways for that group.

- **both**: Jaccard first, then collapsePathways.
