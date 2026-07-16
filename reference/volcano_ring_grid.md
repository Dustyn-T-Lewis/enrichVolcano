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
  subtitles = NULL,
  nrow = NULL,
  ncol = NULL,
  tag_levels = "A",
  guides = "collect",
  panel_spacing = 1.5,
  panel_margin = 2,
  label_headroom = 1.1,
  legend_position = c("bottom", "right", "none"),
  legend_width = 26,
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

- subtitles:

  Optional per-ring subtitles. Either a vector parallel to `contrasts`,
  or a vector named by contrast. `NULL` draws no subtitles.

- nrow, ncol:

  Outer layout dims; forwarded to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html).

- tag_levels:

  Panel-tag scheme; forwarded to
  [`patchwork::plot_annotation()`](https://patchwork.data-imaginist.com/reference/plot_annotation.html).

- guides:

  Patchwork `guides` argument; default `"collect"` collects the shared
  NES legend.

- panel_spacing:

  Gutter between adjacent panels, in millimetres (default 1.5). Applied
  as half on each panel edge, so neighbours sit `panel_spacing` apart.
  Lower it to pack rings closer.

- panel_margin:

  Outer margin around the whole grid, in millimetres (default 2). Trims
  the dead frame around the assembled figure.

- label_headroom:

  Radial room reserved for pathway labels inside each panel (default
  1.1, looser than the single-ring default of 0.5 so wide boxes stay
  enclosed when panels are packed tight). Forwarded to
  [`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md);
  raise it if labels clip, lower it to enlarge the rings.

- legend_position:

  Placement of the collected NES legend: `"bottom"` (default, recovers
  the right-hand gap), `"right"`, or `"none"`.

- legend_width:

  Length of the NES colourbar long axis, in millimetres (default 26,
  tuned for the bottom bar). Sets the key width when the legend is
  horizontal, the key height when vertical; a side legend usually wants
  a larger value (~40).

- ...:

  Forwarded to each
  [`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
  call (e.g. `gene_col`, `padj_col`, `magnitude`, `theme`).

## Value

An S3 object `c("volcano_ring_grid", "list")` with elements `$plot`
(patchwork) and `$data` (list of `list(volc, enrich)` pairs).

## Examples

``` r
da <- read.csv(system.file("extdata", "examples", "yvo_da.csv",
  package = "enrichVolcano"
))
en <- read.csv(system.file("extdata", "examples", "yvo_enrichment.csv",
  package = "enrichVolcano"
))
names(da)[names(da) == "adj.P.Val"] <- "padj"

# up to eight GO-BP terms per contrast, one composite each
en_go <- en[en$database == "GO_BP" & en$padj < 0.01, ]
en_go <- en_go[order(en_go$padj), ]
en_go <- en_go[!duplicated(en_go[c("contrast", "pathway")]), ]
en_go <- do.call(rbind, lapply(split(en_go, en_go$contrast), head, 8))

g <- volcano_ring_grid(da, en_go, contrasts = c("Training_Young", "Training_Old"))
#> No tick-line column found; tick lines off.
#> No tick-line column found; tick lines off.
g$plot
```
