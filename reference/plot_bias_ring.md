# Every significant term of one contrast around its volcano

Draws all significant terms of a contrast as a ring, without names, to
show whether the contrast leans up or down. Every arc gets the same
angle, so the up and down halves grow with their term counts. Fill and
arc height are the score over the strongest drawn score: the strongest
term is full colour and full height, and every other term is a fraction
of it.

## Usage

``` r
plot_bias_ring(
  da,
  enrichment,
  contrast = NULL,
  databases = NULL,
  collapse = TRUE,
  term_threshold = 0.05,
  title = NULL,
  subtitle = NULL,
  theme = plot_theme(),
  ...
)
```

## Arguments

- da:

  DA results from
  [`as_da()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_da.md);
  the rows for `contrast` are drawn. Points are called significant on
  `padj`, or on `p` when `padj` is empty. To colour by another
  statistic, such as a pi-value, pass it to
  [`as_da()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_da.md)
  as `padj`.

- enrichment:

  An enrichment object from
  [`as_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_enrichment.md)
  or
  [`run_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/run_enrichment.md).

- contrast:

  Which contrast of `enrichment` to draw. Needed only when it holds more
  than one.

- databases:

  Collections to draw from, matched against the `database` column, e.g.
  `c("Hallmark", "GO Slim")`. `NULL` (default) draws from all. Skipped,
  with a note, when the object has no database labels.

- collapse:

  Hide terms flagged `"redundant"` by
  [`dedup_terms()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/dedup_terms.md).

- term_threshold:

  Terms need `padj` below this to be drawn.

- title:

  Title.

- subtitle:

  Subtitle. `NULL` reports the up and down term counts.

- theme:

  Output of
  [`plot_theme()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_theme.md).

- ...:

  Other
  [`plot_volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_volcano_ring.md)
  arguments that set the layout, such as `ring_radius`,
  `arc_height_range`, `p_threshold` or `point_size`.

## Value

A ggplot.

## Examples

``` r
da <- read_example("yvo")$da
#> 2103 of 2106 accessions mapped to Homo sapiens symbols.
#> 2106 of 2106 matrix proteins have DA results.
ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
  package = "enrichVolcano"
)))
#> Reading "fgsea" results.
plot_bias_ring(da, ex, contrast = "Aging", title = "Aging")
```
