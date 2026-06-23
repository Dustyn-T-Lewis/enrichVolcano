# Compose a grid of `volcano_ring()` plots, one per contrast

Compose a grid of
[`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
plots, one per contrast

## Usage

``` r
volcano_ring_grid(
  volc_dfs,
  enrich_dfs,
  contrasts = NULL,
  nrow = NULL,
  ncol = NULL,
  tag_levels = "A",
  guides = "collect",
  ...
)
```

## Arguments

- volc_dfs:

  Named list of tidy DA tibbles, one per contrast. A single data.frame
  carrying a `contrast` column is also accepted and is split.

- enrich_dfs:

  Named list of tidy enrichment tibbles, one per contrast. Same split
  convenience as `volc_dfs`.

- contrasts:

  Character vector of contrast names to include and the order in which
  to draw them. Defaults to `names(volc_dfs)`.

- nrow, ncol:

  Outer layout dims; forwarded to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html).

- tag_levels:

  Panel-tag scheme; forwarded to
  [`patchwork::plot_annotation()`](https://patchwork.data-imaginist.com/reference/plot_annotation.html).

- guides:

  Patchwork `guides` argument; default `"collect"` collects the shared
  NES legend.

- ...:

  Forwarded to each
  [`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
  call (e.g. `gene_col`, `padj_col`, `magnitude`, `theme`).

## Value

An S3 object `c("volcano_ring_grid", "list")` with elements `$plot`
(patchwork) and `$data` (list of `list(volc, enrich)` pairs).
