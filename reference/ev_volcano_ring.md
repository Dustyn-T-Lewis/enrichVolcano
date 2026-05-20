# Build a volcano-in-ring composite for one contrast

Build a volcano-in-ring composite for one contrast

## Usage

``` r
ev_volcano_ring(
  volc_df,
  enrich_df,
  title = NULL,
  volcano_radius = 3.5,
  p_threshold = 0.05,
  point_size = 0.9,
  point_alpha = 0.6,
  label_size = 2.6,
  theme = ev_theme()
)
```

## Arguments

- volc_df:

  Tidy long tibble for a single contrast with `gene`, `logFC`,
  `P.Value`, and `pi_eq2` columns.

- enrich_df:

  Enrichment tibble for the same contrast (`pathway`, `padj`, `NES`,
  `size`, `leading_edge`).

- title:

  Plot title (defaults to the contrast name).

- volcano_radius:

  Inner volcano radius.

- p_threshold:

  Significance cutoff applied to `pi_eq2` for up/down calls.

- point_size, point_alpha:

  Volcano point aesthetics.

- label_size:

  Pathway-label text size.

- theme:

  Output of
  [`ev_theme()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_theme.md).

## Value

A ggplot object (class `c("enrichVolcano", ...)`).
