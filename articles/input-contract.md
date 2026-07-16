# Input contract

[`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
takes two tidy frames: a differential-abundance (DA) frame and an
enrichment frame. Every column the function consumes is named via an
argument with a sensible default, so the path of least resistance is to
use limma + fgsea conventional names. Everything else is a rename at the
call site.

``` r

library(enrichVolcano)

volcano_ring(
  volc_df, enrich_df,

  # volcano frame
  gene_col = "gene", # protein/gene identifier
  logfc_col = "logFC", # log fold change
  pval_col = "P.Value", # raw p
  padj_col = "padj", # adjusted p (also used for enrich_df)

  # enrichment frame
  term_col = "pathway", # pathway / term name
  nes_col = "NES", # signed effect size for the term
  size_col = "size", # only required when magnitude = "size"
  genes_col = NULL, # NULL triggers tick auto-detect (Q4)
  genes_sep = NULL, # NULL defers to the auto-detected default
  p_threshold = 0.05,
  logfc_threshold = 0,
  magnitude = c("neg_log_padj", "size")
)
```

## The volcano frame

[`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
needs `gene`, `logFC`, and a p-value column. If a `padj` column is
present it gates the up/down call; if not, the raw `pval` is used.

| Column      | Type                | Required? | Notes                            |
|-------------|---------------------|-----------|----------------------------------|
| `gene_col`  | character           | yes       | one row per gene per contrast    |
| `logfc_col` | numeric             | yes       | log-scale effect                 |
| `pval_col`  | numeric in \[0, 1\] | yes       | raw or moderated                 |
| `padj_col`  | numeric in \[0, 1\] | optional  | when present, gates significance |

### Cross-walk

| Upstream | `gene_col` | `logfc_col` | `pval_col` | `padj_col` | Note |
|----|----|----|----|----|----|
| limma `topTable()` | `Gene` or row name | `logFC` | `P.Value` | `adj.P.Val` | rename `adj.P.Val` to `padj`, or pass `padj_col = "adj.P.Val"` |
| DESeq2 `results()` | row name | `log2FoldChange` | `pvalue` | `padj` | shift row name into a `gene` column |
| edgeR `topTags()` | row name | `logFC` | `PValue` | `FDR` | rename `FDR` to `padj` |
| DEP / proteoDA | `name` | `_ratio` columns | `_p.val` columns | `_p.adj` columns | wide → long first |
| MaxQuant / Perseus | `Gene names` | `t-test difference` | `-LOG(t-test p-value)` | none | back-transform `-log10 p` first |

## The enrichment frame

The enrichment frame is one row per pathway. `NES` is the signed
effect-size column (positive for up-regulated, negative for
down-regulated); `padj` is the multiple-testing-corrected p; `size` is
the gene-set size used by the test.

| Column | Type | Required? | Notes |
|----|----|----|----|
| `term_col` | character | yes | pathway / term name |
| `nes_col` | numeric | yes | signed; NA arcs are dropped |
| `padj_col` | numeric in \[0, 1\] | yes | shared name with the volcano frame |
| `size_col` | numeric | only when `magnitude = "size"` | gene-set size |
| `genes_col` | character or list | optional | leading-edge genes |

### Cross-walk

| Upstream | `term_col` | `nes_col` | `padj_col` | `size_col` | `genes_col` (Q4) |
|----|----|----|----|----|----|
| `fgsea::fgseaMultilevel()` | `pathway` | `NES` | `padj` | `size` | `leadingEdge` (list-col, auto) |
| `clusterProfiler::gseGO()` | `Description` or `ID` | `NES` | `p.adjust` | `setSize` | `core_enrichment` (“/” sep, auto) |
| `enrichR::enrichr()` | `Term` | derived from `Odds.Ratio` | `Adjusted.P.value` | none for ORA | `Genes` (“;” sep, auto) |
| ORA via clusterProfiler | `Description` | derive from over/under | `p.adjust` | `Count` | `geneID` (“/” sep) |

Pass column overrides when the names don’t match — for example with
clusterProfiler GSEA output:

``` r

volcano_ring(
  da, gseaResult@result,
  term_col = "Description",
  padj_col = "p.adjust",
  size_col = "setSize",
  genes_col = "core_enrichment",
  genes_sep = "/"
)
```

## Tick-column auto-detect (Q4)

When `genes_col = NULL` (the default),
[`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
walks four conventional column names in order:

1.  `leading_edge` (`;`-separated string)
2.  `leadingEdge` (list-column, fgsea native)
3.  `core_enrichment` (`/`-separated string, clusterProfiler GSEA)
4.  `Genes` (`;`-separated string, enrichR)

First match wins. A `cli_alert_info()` names the column. No match
silently turns tick lines off — the rings still draw, just without the
radial tick segments inside the volcano.

Override with `genes_col = "my_column"` and, for a string column, a
matching `genes_sep` (defaults to `";"` when the column is a string).

## Renaming strategy

Three patterns work, in increasing order of intrusiveness:

1.  **Argument-only renames**:
    `volcano_ring(..., padj_col = "adj.P.Val")`. No frame mutation, good
    for one-off plots.
2.  **Boundary rename**:
    `names(da)[names(da) == "adj.P.Val"] <- "padj"`. Useful when several
    functions in the pipeline expect the same name.
3.  **Built-in column shim**: write a tiny adapter `to_evolc(df)` that
    returns the renamed frame. Worth doing when you have a
    project-specific upstream format you reuse a lot.

When the DA and enrichment sides disagree on the `padj` column name
(e.g. `adj.P.Val` vs `padj`), one of them has to be renamed at the
boundary — `padj_col` is a single shared argument by design.

## Worked example: DEP + fgsea

`DEP::get_results()` returns a wide tibble with one row per protein and
one group of columns per contrast (`<contrast>_ratio`,
`<contrast>_p.val`, `<contrast>_p.adj`).
[`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
wants one row per protein **per contrast**, so a single
`tidyr::pivot_longer()` reshapes it:

``` r

library(tidyr)
library(dplyr)

da <- DEP::get_results(dep) |>
  select(name, ends_with(c("_ratio", "_p.val", "_p.adj"))) |>
  pivot_longer(
    -name,
    names_to = c("contrast", ".value"),
    names_pattern = "(.*)_(ratio|p\\.val|p\\.adj)"
  ) |>
  rename(gene = name, logFC = ratio, P.Value = p.val, padj = p.adj)
```

Enrichment from `fgsea` is already in the shape
[`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
expects — no reshape, no rename:

``` r

ranks <- da |>
  filter(contrast == "treated_vs_control") |>
  with(setNames(logFC, gene))

en <- fgsea::fgseaMultilevel(pathways = msigdbr_list, stats = ranks) |>
  mutate(contrast = "treated_vs_control")

volcano_ring(
  volc_df   = filter(da, contrast == "treated_vs_control"),
  enrich_df = en
)
```

## NA + edge-case policy

| Case                  | Behaviour                          |
|-----------------------|------------------------------------|
| logFC or P.Value NA   | row dropped                        |
| every NES NA          | error (cannot colour arcs)         |
| pathway NA            | row dropped from the ring          |
| padj NA in enrichment | filled with 1                      |
| NES NA in enrichment  | filled with 0                      |
| duplicate term        | warning + dedup by lowest padj     |
| no tick column found  | tick lines silently off, one alert |
| p outside \[0, 1\]    | error                              |

## Upstream tools we consume

The package is a thin plot wrapper around tidy outputs from existing
enrichment tooling:

- fgsea — Korotkevich G et al. 2021, DOI 10.1101/060012
- clusterProfiler — Wu T et al. 2021 *Innovation*
- enrichR — Chen EY et al. 2013 *BMC Bioinformatics*

For methodology when *choosing* which enrichment to run, see:

- Xiao Y et al. 2014 *Bioinformatics* 30(6):801-7, PMID 22321699
  (π-value)
- Timmons JA et al. 2015 *Genome Biol* 16:186, PMID 26346307 (ORA
  universe)
- Reimand J et al. 2019 *Nat Protoc* 14:482-517, PMID 30664679
  (protocol)
- Wijesooriya K et al. 2022 *PLoS Comput Biol* 18:e1009935, PMID
  35263338 (background bias)
