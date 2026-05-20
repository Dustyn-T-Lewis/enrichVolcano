# Build an enrichment ring ggplot for one contrast

Uses
[`ggforce::geom_arc_bar`](https://ggforce.data-imaginist.com/reference/geom_arc_bar.html)
to display the top-N pathways as arc segments. Magnitude maps to arc
thickness, color to NES (or other), terms split by up/down NES
direction.

## Usage

``` r
ring_plot(
  enrich_result,
  contrast,
  max_terms = 10,
  order_by = c("padj", "size", "NES", "alphabetic"),
  order_dir = c("asc", "desc"),
  magnitude = c("neg_log_padj", "nes", "gene_count"),
  color = c("nes", "neg_log_padj", "direction"),
  split_by_direction = TRUE,
  label_clean_fn = NULL,
  theme = ev_theme()
)
```

## Arguments

- enrich_result:

  Tibble from
  [`ev_collapse()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_collapse.md)
  or
  [`ev_enrich()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_enrich.md).

- contrast:

  Single contrast name.

- max_terms:

  Cap on terms shown (default 10; warn at 12, abort at 25).

- order_by:

  `"padj"` (default), `"size"`, `"NES"`, or `"alphabetic"`.

- order_dir:

  `"asc"` or `"desc"`.

- magnitude:

  What arc-thickness encodes: `"neg_log_padj"` (default), `"nes"`, or
  `"gene_count"`.

- color:

  Arc color encoding: `"nes"` (default), `"neg_log_padj"`, or
  `"direction"`.

- split_by_direction:

  Split up- and down-NES into two arc groups.

- label_clean_fn:

  Optional user fn applied to pathway names.

- theme:

  Output of
  [`ev_theme()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_theme.md).

## Value

ggplot object.
