# Fonts and colours for the plots

Sets the font and the colours that
[`plot_volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_volcano_ring.md)
and
[`plot_scatter()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_scatter.md)
use: the up, down and non-significant point colours and the diverging
ramp that fills arcs by score.

## Usage

``` r
plot_theme(
  base_size = 11,
  base_family = "",
  palette = c("default", "viridis", "okabe"),
  up = NULL,
  down = NULL,
  ns = NULL,
  score_colours = NULL,
  score_limits = NULL,
  score_stops = NULL
)
```

## Arguments

- base_size:

  Base font size.

- base_family:

  Base font family. `""` uses the ggplot2 default.

- palette:

  `"default"` (red and blue), `"viridis"` (magma) or `"okabe"`
  (Okabe-Ito, safe for colour-blind readers). `up`, `down`, `ns` and
  `score_colours` override single parts of it.

- up, down, ns:

  Colours of up-regulated, down-regulated and non-significant points.

- score_colours:

  Colours of the score ramp, low to high. Without `score_stops` they
  spread evenly across `score_limits`, or `c(-3, 3)`.

- score_limits:

  Length-2 limits of the score scale. `NULL` lets each plot choose; see
  `score_limits` in
  [`plot_volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_volcano_ring.md).

- score_stops:

  Score values at which each ramp colour sits, one per colour.

## Value

A list with `base_size`, `base_family`, `palette` and `score_limits`,
passed to the plots as `theme`.

## Examples

``` r
th <- plot_theme(base_size = 11)
th$palette$up
#> [1] "#D6604D"

th <- plot_theme(up = "#B2182B", down = "#2166AC", ns = "grey80")
```
