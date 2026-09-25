# Compare two contrasts term by term

Plots each term's score in one contrast against its score in another.
Terms in the top-right and bottom-left quadrants move the same way in
both contrasts; terms in the other two move in opposite directions.

## Usage

``` r
plot_scatter(
  enrichment,
  x,
  y,
  comparison = c("concordance", "reversal"),
  databases = NULL,
  collapse = TRUE,
  term_threshold = 0.05,
  colour_by = "significance",
  shape_by = "database",
  label_min_size = 15,
  label_n = 20,
  labels = NULL,
  theme = plot_theme()
)
```

## Arguments

- enrichment:

  An enrichment object from
  [`as_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_enrichment.md)
  or
  [`run_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/run_enrichment.md)
  that holds both contrasts.

- x, y:

  Contrast names for the horizontal and vertical axes.

- comparison:

  `"concordance"` or `"reversal"`; see Description.

- databases:

  Collections to draw from, e.g. `c("Hallmark", "GO Slim")`; `NULL`
  (default) draws from all. Skipped, with a note, when the object has no
  database labels.

- collapse:

  Hide terms that
  [`dedup_terms()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/dedup_terms.md)
  flagged redundant in one contrast and kept in neither.

- term_threshold:

  A term is significant in a contrast when its `padj` is below this.

- colour_by:

  Column that fills significant points: `"significance"` (significant in
  `x` only, `y` only, or both), any per-term column such as
  `"database"`, or `NULL` for one colour. Per-term columns are read from
  the `x` contrast.

- shape_by:

  Column mapped to point shape (up to five values), or `NULL`.

- label_min_size:

  Smallest gene set that gets a label.

- label_n:

  Most term labels drawn.

- labels:

  Your own display names, as in
  [`plot_volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_volcano_ring.md).

- theme:

  Output of
  [`plot_theme()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_theme.md);
  supplies the base font.

## Value

A ggplot.

## Details

`comparison` sets the question the figure asks:

- `"concordance"`: do two conditions regulate the same pathways, e.g.
  training in young against training in old? Same-direction quadrants
  are "concordant", the reference line is \\y = x\\, and the subtitle
  reports the share of significant terms that are concordant.

- `"reversal"`: does one condition undo another, e.g. aging against
  training? Opposite-direction quadrants are "reversed", same-direction
  ones "exacerbated", the reference line is \\y = -x\\, and the subtitle
  reports the share reversed.

## What is drawn

- Points: every term scored in both contrasts. Non-significant terms are
  small and grey; significant ones are sized by gene-set size and filled
  by `colour_by`.

- Shading: same-direction quadrants are tinted, the others left white.

- Corner labels: each quadrant's name and its count of significant
  terms.

- Dashed line: \\y = x\\, or \\y = -x\\ for a reversal.

- Labels: the most significant terms with at least `label_min_size`
  genes.

- Subtitle: Spearman's \\\rho\\ across all plotted terms, with a 95%
  Fieller interval when the correlation package is installed, its
  p-value, the share of significant terms that are concordant (or
  reversed), and the term counts.

A term is significant in a contrast when its `padj` is below
`term_threshold`. If the object has no adjusted p-values at all, nominal
p-values are used and a note says so.

## References

Fieller EC, Hartley HO, Pearson ES (1957). Tests for rank correlation
coefficients. I. Biometrika 44(3-4):470-481.
[doi:10.1093/biomet/44.3-4.470](https://doi.org/10.1093/biomet/44.3-4.470)

Makowski D, Ben-Shachar MS, Patil I, Luedecke D (2020). Methods and
algorithms for correlation analysis in R. Journal of Open Source
Software 5(51):2306.
[doi:10.21105/joss.02306](https://doi.org/10.21105/joss.02306)

## Examples

``` r
ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
  package = "enrichVolcano"
)))
#> Reading "fgsea" results.
plot_scatter(ex, "Training_Young", "Training_Old")
#> 2 terms appear in only one contrast and are left out.

plot_scatter(ex, "Aging", "Training_Old", comparison = "reversal")
#> 2 terms appear in only one contrast and are left out.
```
