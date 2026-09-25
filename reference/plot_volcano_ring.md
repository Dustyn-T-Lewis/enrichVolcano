# Volcano-in-ring composite for one contrast

Draws a differential-abundance volcano embedded in a ring of enrichment
terms. Score-coloured arcs sit around the volcano; tick lines drop from
the arcs to each term's leading-edge genes inside the volcano.

## Usage

``` r
plot_volcano_ring(
  da,
  enrichment,
  contrast = NULL,
  databases = NULL,
  collapse = TRUE,
  term_threshold = 0.05,
  n_terms = 12,
  terms = NULL,
  labels = NULL,
  p_threshold = 0.05,
  logfc_threshold = 0,
  title = NULL,
  subtitle = NULL,
  tag = NULL,
  volcano_radius = 4,
  x_scale = 1,
  y_scale = 1,
  ring_radius = 4.8,
  ring_thickness = 0.55,
  tick_width = 0.3,
  label_headroom = 0.5,
  disc_colour = NULL,
  score_limits = NULL,
  magnitude = NULL,
  arc_order = c("padj", "score"),
  arc_height_range = c(0.4, 1.6),
  show_counts = TRUE,
  point_size = 1.1,
  point_alpha = 0.85,
  label_size = 2.8,
  label_gap = 0.6,
  count_size = 2.4,
  count_x_mult = 0.7,
  count_y_mult = 0.7,
  axis_size = 2.2,
  label_mode = c("none", "top_per_direction", "by_significance", "by_genes"),
  label_n = 5,
  label_rank_by = c("significance", "logfc"),
  label_genes = NULL,
  theme = plot_theme()
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

- n_terms:

  Most unique terms drawn, counted across both directions.

- terms:

  Optional character vector of exact term names to draw instead.

- labels:

  Your own display names, as a character vector named by term, such as
  `c(HALLMARK_OXIDATIVE_PHOSPHORYLATION = "OXPHOS")`. Terms not named
  keep
  [`clean_label()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/clean_label.md).
  A name without `\n` is wrapped like the others.

- p_threshold:

  Significance cutoff for volcano points.

- logfc_threshold:

  Effect-size cutoff; a point is called up/down only when
  `abs(logFC) >= logfc_threshold` as well as significant.

- title, subtitle, tag:

  Title, subtitle and panel tag.

- volcano_radius:

  Radius of the volcano in plot units.

- x_scale:

  Horizontal compression of the point cloud (default 1). Values below 1
  pull points toward the fold-change axis so the widest points clear the
  enrichment ring; the up/down axis annotations are unaffected.

- y_scale:

  Vertical compression of the point cloud (default 1), anchored at the
  fold-change axis. Values below 1 lower the tallest points so they
  clear the `-log10 p` label at the top of the volcano.

- ring_radius:

  Inner radius of the enrichment ring (default 4.8), where the
  leading-edge tick band begins. Raise it to widen the central breathing
  gap around the volcano, lower it to close it. Keep it above
  `volcano_radius * 0.92` so the point cloud clears the ring.

- ring_thickness:

  Radial width of the tick band between `ring_radius` and the foot of
  the coloured arcs (default 0.55). This is the length of the
  leading-edge ticks; widen it to make ticks easier to read.

- tick_width:

  Line width of the leading-edge ticks (default 0.3).

- label_headroom:

  Extra radial room (data units, default 0.5) reserved beyond the
  outermost pathway label so its box stays enclosed within the square
  panel rather than clipping or spilling into a neighbour. Raise it when
  wide label boxes are clipped; lower it to pack the ring tighter.

- disc_colour:

  Optional fill for a tinted central disc.

- score_limits:

  Length-2 numeric or `NULL`; limits of the arc fill scale. `NULL` uses
  `c(-3, 3)` for NES, and otherwise spans the largest absolute
  significant score in the whole object, so every contrast shares a
  scale.

- magnitude:

  What arc height encodes: `"neg_log_padj"` or `"size"`. `NULL` picks
  `"neg_log_padj"` for NES and `"size"` otherwise, so fill and height
  never repeat the same number for fry or camera results; without set
  sizes it falls back to `"neg_log_padj"`.

- arc_order:

  Angular order of arcs within each up/down half: `"padj"` (default,
  lowest FDR first) or `"score"` (strongest absolute score first). The
  up/down split itself is always by direction.

- arc_height_range:

  Length-2 numeric `c(min, max)` for the shortest and tallest arc; widen
  it to exaggerate the magnitude encoding.

- show_counts:

  Draw the up/down significant-point count badges.

- point_size, point_alpha:

  Volcano point size (default 1.1) and opacity.

- label_size:

  Text size of the term labels and the point labels.

- label_gap:

  Radial gap between each arc's outer top and its own label (default
  0.6). Anchoring per-arc keeps every leader line the same short length
  regardless of arc height; the label box grows outward from this point
  so it never overlaps the arc. Widen it to lengthen all leaders.

- count_size:

  Text size of the up/down count badges (default 2.4).

- count_x_mult, count_y_mult:

  Badge position as a fraction of the volcano radius (default 0.7).

- axis_size:

  Text size of the `up`/`down`/`log2 FC`/`-log10 p` axis annotations
  (default 2.2).

- label_mode:

  Which points get a gene label: `"none"`, `"top_per_direction"`
  (`label_n` up and `label_n` down), `"by_significance"` (`label_n` in
  total) or `"by_genes"` (the proteins in `label_genes`).

- label_n:

  How many points the two top modes label.

- label_rank_by:

  Rank points for the top modes by `"significance"` (smallest `p`) or
  `"logfc"` (largest absolute fold change).

- label_genes:

  Gene symbols or accessions to label with `"by_genes"`.

- theme:

  Output of
  [`plot_theme()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_theme.md).

## Value

A ggplot.

## Which terms ring the volcano

From the chosen contrast, terms in `databases` are kept, redundant ones
are hidden when `collapse = TRUE` (see
[`dedup_terms()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/dedup_terms.md)),
and the `n_terms` terms with the smallest `padj` below `term_threshold`
are drawn, whatever their direction, so the up and down halves need not
match. A term found in two collections is drawn once. `terms` overrides
all of this with a hand-picked set.

## Examples

``` r
da <- read_example("yvo")$da
#> 2103 of 2106 accessions mapped to Homo sapiens symbols.
#> 2106 of 2106 matrix proteins have DA results.
ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
  package = "enrichVolcano"
)))
#> Reading "fgsea" results.

plot_volcano_ring(da, ex, contrast = "Training_Young", title = "Training_Young")
```
