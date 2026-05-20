# Introduction to enrichVolcano

## What this package does

`enrichVolcano` builds one figure with two panels: a volcano plot of
differential abundance results and a ring plot of pathway enrichment for
the same contrast. Input can come from `proteoDA`, `DEP`, `limma`,
`DESeq2`, `edgeR`, `MSstats`, `proDA`, `DEqMS`, MaxQuant, or Perseus;
[`ev_validate()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_validate.md)
auto-detects the format and returns a tidy long tibble. The hero
function
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
returns a `patchwork` object that you can save with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
or further compose with other plots.

``` r

library(enrichVolcano)
```

## A 30-second example

Load the bundled YvO (young vs. old skeletal muscle) result set. It
ships with the package after install; during development you can read it
from the test fixtures.

``` r

fixture <- system.file(
  "extdata", "examples", "yvo_tidy.rds",
  package = "enrichVolcano"
)
if (!nzchar(fixture)) {
  fixture <- "../tests/testthat/fixtures/yvo_tidy.rds"
}
yvo <- readRDS(fixture)
yvo <- yvo[!is.na(yvo$gene) & !is.na(yvo$contrast), ]
head(yvo, 3)
#>    gene contrast     logFC         t          B      P.Value   adj.P.Val
#> 1 LMOD2    Aging -2.347025 -4.532596  2.1657602 3.645698e-05 0.001783001
#> 2 MYBPH    Aging -2.271899 -4.048600  0.7482425 1.685745e-04 0.004726829
#> 3 ACTN3    Aging -2.120998 -3.211231 -1.6316218 2.227110e-03 0.024142332
#>       pi_score sig_pi
#> 1 3.831672e-11     -1
#> 2 2.677030e-09     -1
#> 3 2.369003e-06     -1
unique(yvo$contrast)
#> [1] "Aging"          "Interaction"    "Training_Old"   "Training_Young"
```

The fixture has 400 protein rows across four contrasts (Aging,
Training_Young, Training_Old, Interaction). Each row already has
`logFC`, `P.Value`, and `adj.P.Val`.

The default registered databases (`hallmark`, `reactome`, `go_bp`) fetch
from `msigdbr`, which needs an internet connection. For a self-contained
example we build a small mock pathway list from the genes in the
fixture.

``` r

set.seed(42)
genes <- unique(yvo$gene)
mock_msigdb <- list(
  HALLMARK_GLYCOLYSIS         = sample(genes, 30),
  HALLMARK_OXIDATIVE_PHOS     = sample(genes, 25),
  HALLMARK_APOPTOSIS          = sample(genes, 20),
  HALLMARK_MYOGENESIS         = sample(genes, 35),
  REACTOME_CITRIC_ACID_CYCLE  = sample(genes, 15)
)
```

Now call
[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
with one contrast and the mock database:

``` r

fig <- enrich_volcano(
  data       = yvo,
  contrast   = "Aging",
  databases  = list(mock = mock_msigdb),
  enrich_padj = 1,
  enrich_mode = "fgsea",
  ring   = list(max_terms = 5)
)
fig
#> 
#> ── enrichVolcano ──
#> 
#> Call: `enrich_volcano(data = yvo, contrast = "Aging", databases = list(mock =
#> mock_msigdb), enrich_mode = "fgsea", enrich_padj = 1, ring = list(max_terms =
#> 5))`
#> Contrasts: 4
#> Proteins: 400
#> Significant pathways (post-dedup): 0
```

![Volcano panel above an enrichment ring panel for the Aging
contrast.](enrichVolcano_files/figure-html/hero-call-1.png)

`enrich_padj = 1` keeps every pathway in the ring for the demo; in
practice you would use the default 0.05. The `mock_msigdb` here stands
in for the real Hallmark or Reactome collections.

## Anatomy of the output

The volcano sits at the centre: `logFC` on the x-axis, `-log10(pi_eq2)`
on the y-axis. A point is called up or down when its pi-score clears
`p_threshold` (default 0.05) and its `logFC` is positive or negative;
the panel reports the up and down counts in coloured boxes rather than
labelling individual genes. (For a labelled standalone volcano with
`ggrepel`, call
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md)
— see
[`vignette("customising")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/customising.md).)

The enrichment ring wraps the volcano, drawn with
[`ggforce::geom_arc_bar`](https://ggforce.data-imaginist.com/reference/geom_arc_bar.html)
in a fixed Cartesian layout (no `coord_polar`). Up-regulated pathways
cluster on the right (centred at 90°) and down-regulated on the left
(270°), so direction reads at a glance. Each arc’s outer radius encodes
`-log10(padj)` (taller = more significant) and its fill encodes NES on a
colourblind-safe blue-white-red diverging scale. Short radial ticks
under each arc mark its leading-edge genes, coloured by the sign of
their fold change.

The contrast name is the title. The composite is a `patchwork` object,
so it supports the `+`, `/`, and `|` operators for further layout.

Two attributes hold the intermediate state:

``` r

str(attr(fig, "ev_data"), max.level = 1)
#> List of 4
#>  $ validated_input: tibble [400 × 9] (S3: tbl_df/tbl/data.frame)
#>   ..- attr(*, "ev_source")= chr "limma"
#>  $ pi_scores      : tibble [400 × 10] (S3: tbl_df/tbl/data.frame)
#>   ..- attr(*, "ev_source")= chr "limma"
#>  $ enrichment     : tibble [0 × 13] (S3: tbl_df/tbl/data.frame)
#>   ..- attr(*, "ev_pathways")=List of 1
#>   ..- attr(*, "ev_stats")=List of 1
#>  $ dedup_result   : tibble [0 × 13] (S3: tbl_df/tbl/data.frame)
#>   ..- attr(*, "ev_pathways")=List of 1
#>   ..- attr(*, "ev_stats")=List of 1
attr(fig, "ev_call")
#> enrich_volcano(data = yvo, contrast = "Aging", databases = list(mock = mock_msigdb), 
#>     enrich_mode = "fgsea", enrich_padj = 1, ring = list(max_terms = 5))
```

`ev_data` is a list with four tibbles: the validated input, the pi-score
table, the raw
[`ev_enrich()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_enrich.md)
output, and the post-dedup enrichment table. `ev_call` is the literal
[`match.call()`](https://rdrr.io/r/base/match.call.html), so a reviewer
can re-run the exact figure from the object alone.

## The pipeline under the hood

[`enrich_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrich_volcano.md)
chains the exported pipeline stages and then draws the canonical
composite:

1.  `ev_validate(data)` detects the input format and renames columns to
    `gene`, `contrast`, `logFC`, `P.Value`, `adj.P.Val`.
2.  `pi_score(data, variant = "eq2")` computes the Xiao 2014 pi-score
    and writes it to `pi_eq2`.
3.  `adjust_p(data, method = "BH")` recomputes `adj.P.Val` from the raw
    p-values within each contrast.
4.  `ev_enrich(data, contrast, databases)` runs `fgsea` (and optionally
    `fora`) against the requested pathway sets, ranking genes by the
    signed statistic `sign(logFC) * -log10(P)`.
5.  `ev_collapse(enrich, method = "jaccard")` removes redundant pathways
    with overlapping leading-edge gene sets.
6.  `ev_volcano_ring(volcano_sub, enrich_sub)` draws one volcano-in-ring
    composite per contrast;
    [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
    arranges multiple contrasts and attaches `ev_data` and `ev_call`.

Stages 1-5 are exported, so a real analysis can stop between them — a
common pattern is to compute pi-scores once, write the table to disk,
and re-enter at
[`ev_enrich()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_enrich.md)
with different database choices. For a panel layout you assemble
yourself, the package also ships the standalone trio
[`ev_volcano()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_volcano.md) +
[`ring_plot()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ring_plot.md) +
[`ev_compose()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_compose.md),
which expose the tunable label and ring-encoding options the opinionated
composite keeps fixed (see
[`vignette("customising")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/customising.md)).

## What next

- [`vignette("scoring")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/scoring.md)
  covers the pi-score variants and the four p-value adjustment methods.
- [`vignette("pathway-dedup")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/pathway-dedup.md)
  explains Jaccard collapsing, scope (within vs. across databases), and
  when to use
  [`fgsea::collapsePathways`](https://rdrr.io/pkg/fgsea/man/collapsePathways.html)
  instead.
- [`vignette("databases")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/databases.md)
  lists the registered gene-set sources and shows how to add a custom
  GMT.
- [`vignette("customising")`](https://Dustyn-T-Lewis.github.io/enrichVolcano/articles/customising.md)
  walks through theme overrides, palette swaps, and per-panel ggplot
  edits.
- `vignette("faq")` answers questions about MaxQuant ratios, contrast
  naming, and large figure layouts.
