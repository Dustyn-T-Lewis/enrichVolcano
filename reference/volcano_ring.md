# Volcano-in-ring composite for one contrast

Draws a differential-abundance volcano embedded in a ring of enrichment
terms. NES-coloured arcs sit around the volcano; tick lines drop from
the arcs to each pathway's leading-edge genes inside the volcano.

## Usage

``` r
volcano_ring(
  volc_df,
  enrich_df,
  gene_col = "gene",
  logfc_col = "logFC",
  pval_col = "P.Value",
  padj_col = "padj",
  term_col = "pathway",
  nes_col = "NES",
  size_col = "size",
  genes_col = NULL,
  genes_sep = NULL,
  p_threshold = 0.05,
  logfc_threshold = 0,
  title = NULL,
  subtitle = NULL,
  tag = NULL,
  volcano_radius = 3.5,
  disc_color = NULL,
  nes_limits = NULL,
  magnitude = c("neg_log_padj", "size"),
  point_size = 1.6,
  point_alpha = 0.85,
  label_size = 2.8,
  count_x_mult = 0.9,
  count_y_mult = 0.9,
  label_mode = c("none", "top_per_direction", "by_significance", "by_genes"),
  label_n = 5,
  label_rank_by = c("significance", "logfc"),
  label_genes = NULL,
  theme = volcano_ring_theme()
)
```

## Arguments

- volc_df:

  Tidy DA tibble for a single contrast.

- enrich_df:

  Tidy enrichment tibble for the same contrast.

- gene_col, logfc_col, pval_col, padj_col:

  Column names in `volc_df`.

- term_col, nes_col, size_col:

  Column names in `enrich_df`. `padj_col` is reused for the
  enrichment-side adjusted-p column.

- genes_col:

  Optional column name in `enrich_df` carrying leading-edge genes.
  `NULL` triggers auto-detect (Q4).

- genes_sep:

  Separator for the leading-edge string. `NULL` defers to the
  auto-detected default.

- p_threshold:

  Cutoff applied to `padj_col` in `volc_df`.

- logfc_threshold:

  Effect-size cutoff; a point is called up/down only when
  `abs(logFC) >= logfc_threshold` as well as significant.

- title, subtitle, tag:

  Plot text.

- volcano_radius:

  Inner volcano radius.

- disc_color:

  Optional fill for a tinted central disc.

- nes_limits:

  Length-2 numeric or `NULL`; defaults to `c(-3, 3)`.

- magnitude:

  `"neg_log_padj"` (default) or `"size"`; controls arc thickness
  encoding.

- point_size, point_alpha, label_size, count_x_mult, count_y_mult:

  Aesthetic controls for the volcano interior.

- label_mode, label_n, label_rank_by, label_genes:

  Volcano point labels.

- theme:

  Output of
  [`volcano_ring_theme()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring_theme.md).

## Value

A ggplot.

## Details

Column-naming arguments default to limma + fgsea conventions. The tick
column is auto-detected from `leading_edge` / `leadingEdge` /
`core_enrichment` / `Genes` unless `genes_col` is supplied.
