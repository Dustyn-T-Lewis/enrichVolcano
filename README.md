
# enrichVolcano <a href="https://Dustyn-T-Lewis.github.io/enrichVolcano/"><img src="man/figures/logo.png" align="right" height="139" alt="enrichVolcano website" /></a>

> Figures from enrichment results you have already computed: fgsea,
> clusterProfiler, limma’s fry and camera, over-representation, or your
> own.

<img src="man/figures/README-hero-1.png" alt="A volcano plot ringed by enrichment arcs coloured by normalised enrichment score"  />

enrichVolcano draws; it does not run differential abundance or
enrichment. Every workflow has three steps:

``` r
library(enrichVolcano)

x <- as_enrichment(list(Aging = fgsea_result), database = "Hallmark")
x <- dedup(x, gene_sets)
volcano_ring(da_aging, x)
```

1.  **`as_enrichment()`** converts results from any supported tool into
    one validated `enrichment` object and records which test produced
    them.
2.  **`dedup()`** (optional) flags redundant terms so a figure shows
    each signal once. p-values are never changed.
3.  **Plot**: `volcano_ring()` for one contrast, `volcano_ring_grid()`
    for several, `nes_scatter()` to compare two contrasts term by term
    (concordance or reversal).

## Install

``` r
install.packages("remotes")
remotes::install_github("Dustyn-T-Lewis/enrichVolcano")
```

## What goes in

| Tool | Recognised by | Score drawn |
|----|----|----|
| fgsea | its columns | NES |
| clusterProfiler `gseaResult` | its class | NES |
| limma `fry`, `mroast` | their columns | signed −log10(FDR) |
| limma `camera`, `cameraPR` | `enrichment_test =` | signed −log10(FDR) |
| over-representation (enrichR, `enrichGO`, …) | `enrichment_test = "ora"` + a direction column | signed −log10(padj) |
| anything else | `enrichment_test = "custom"` + column names | yours |

fry, mroast and camera report no effect size or leading edge, so their
score is the signed −log10(FDR) and the legend says so.

## What the ring shows

- **Centre:** the volcano, `logFC` against −log10(p).
- **Ring:** one arc per term, up-regulated on top and down-regulated
  below. Fill is the score; height is −log10(padj) for NES and set size
  otherwise.
- **Ticks:** spokes from each arc to its leading-edge genes in the
  volcano.

By default the ring draws up to eight significant terms per direction
from Hallmark and GO slim; `databases`, `n_terms`, `term_threshold` and
`terms` change that. `vignette("enrichVolcano")` walks every step on the
bundled example.

## Cite the methods you used

enrichVolcano only draws; the credit belongs to the analysis:

- fgsea: Korotkevich G et al., bioRxiv, <doi:10.1101/060012>
- clusterProfiler: Wu T et al. 2021, *The Innovation* 2:100141
- camera: Wu D, Smyth GK 2012, *Nucleic Acids Res* 40(17),
  <doi:10.1093/nar/gks461>
- fry and mroast (ROAST): Wu D et al. 2010, *Bioinformatics* 26(17),
  <doi:10.1093/bioinformatics/btq401>
- EnrichmentMap (`dedup()`): Merico D et al. 2010, *PLoS ONE*
  5(11):e13984
- Enrichment protocol: Reimand J et al. 2019, *Nat Protoc* 14:482-517

## Versioning

1.0.0 replaced the column-name arguments with the `enrichment` object;
see `NEWS.md`. Earlier releases stay installable from their tags.

## License

MIT. See `LICENSE`.
