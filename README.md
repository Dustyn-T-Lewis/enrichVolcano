
<!-- README.md is generated from README.Rmd. Please edit that file -->

# enrichVolcano

<!-- badges: start -->

[![R-CMD-check](https://github.com/Dustyn-T-Lewis/enrichVolcano/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Dustyn-T-Lewis/enrichVolcano/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

enrichVolcano builds composite volcano + enrichment ring plots from
differential abundance results. The volcano sits at the centre and the
enrichment terms wrap around it as NES-coloured arcs, one figure per
contrast. Returns a ggplot-compatible object that `ggsave()` saves
directly.

The package accepts twelve common input formats (tidy, wide-suffix, DEP,
proteoDA, limma, DESeq2, edgeR, MSstats, proDA, DEqMS, MaxQuant,
Perseus) via `ev_validate()`, supports two pi-score variants (Xiao 2014
Eq.1 and Eq.2), four p-adjustment methods (BH, Bonferroni, q-value,
IHW), and 19 registered gene-set databases across seven species.

## Installation

``` r
# Development version
remotes::install_github("Dustyn-T-Lewis/enrichVolcano")
```

## A 30-second example

``` r
library(enrichVolcano)

# Load a bundled YvO subsample
data <- readRDS(system.file("extdata", "examples", "yvo_tidy.rds",
                            package = "enrichVolcano"))

# Mock pathway sets (skip msigdbr in this demo)
paths <- list(
  UP_GENES = head(data$gene[data$logFC > 0], 30),
  DN_GENES = head(data$gene[data$logFC < 0], 30)
)

p <- enrich_volcano(
  data, contrast = "Aging",
  databases = list(test = paths),
  p_method = "pi_eq2", p_adjust = "BH",
  enrich_padj = 0.5
)

p
```

The output prints the reconstructible call plus a one-line data summary,
then renders the patchwork. Both panels share the contrast title.

## Why a composed plot

Volcano panels show per-protein effects. Enrichment plots show
per-pathway aggregates. Drawing them separately loses the connection
between a hit and the pathway it drives. The composite preserves the
hit-to-pathway lineage and produces one figure per contrast instead of
two.

## Reproducibility

Every output carries two attributes:

- `attr(p, "ev_call")` — `match.call()` for exact reconstruction
- `attr(p, "ev_data")` — list of intermediate tibbles: validated input,
  pi-scores, enrichment, dedup result

These let you rebuild any panel without re-running the pipeline.

## Pipeline functions

`enrich_volcano()` runs the pipeline below and draws the composite with
`ev_volcano_ring()`. Every stage is exported, so you can run or swap any
one yourself:

``` r
ev_validate()      # input adapter (12 sources)
pi_score()         # Xiao 2014 Eq.1 or Eq.2
adjust_p()         # BH / Bonferroni / qvalue / IHW
ev_enrich()        # fgsea + ORA
ev_collapse()      # Jaccard / collapsePathways dedup
ev_volcano_ring()  # the volcano-in-ring composite
```

For building panels by hand there are also `ev_volcano()` (volcano
only), `ring_plot()` (ring only), and `ev_compose()` (patchwork stitch),
plus `list_databases()` / `database_info()` for the registry and
`ev_theme()` for styling.

## Learn more

- `vignette("enrichVolcano")` — 30-second example + pipeline anatomy
- `vignette("scoring")` — pi-score variants, p-value adjustment methods
- `vignette("pathway-dedup")` — Jaccard vs collapsePathways, scope
  choices
- `vignette("databases")` — 19 databases compared, decision flowchart
- `vignette("customising")` — themes, palettes, multi-contrast layouts
- `vignette("enrichVolcano-faq")` — common questions

## Citation

``` r
citation("enrichVolcano")
```

The package implements the pi-score from Xiao Y et al. (2014)
*Bioinformatics* 30(6):801, fgsea from Korotkevich G et al. (2021)
bioRxiv 060012, and Benjamini-Hochberg FDR from Benjamini & Hochberg
(1995) JRSS B 57(1):289.

## License

MIT
