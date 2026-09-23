# Enrichment results in one validated object

`enrichment` is the single input every plot in the package reads. Build
it with
[`as_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_enrichment.md)
rather than by hand; construct it directly only when you already hold a
table in exactly this shape.

## Usage

``` r
enrichment(results = data.frame(), metadata = list())
```

## Arguments

- results:

  A data frame with one row per term per contrast:

  - `contrast`, `database`, `term`: character. `database` may be `NA`.

  - `score`: numeric. NES for fgsea and clusterProfiler; signed
    \\-\log\_{10}\\(FDR) for limma rotation/competitive tests and ORA.

  - `p`, `padj`: numeric in \[0, 1\]; `NA` allowed.

  - `size`: numeric set size.

  - `direction`: `"up"` or `"down"`, agreeing with the sign of `score`.

  - `leading_edge`: list of character vectors, empty when the test has
    none.

  Extra columns (for example `dedup_status` or your own groupings) are
  kept and can be mapped by the plots.

- metadata:

  A list with `enrichment_test` (the method that produced the results),
  `score_type` (the axis and legend label for `score`), and `dedup`
  (`NULL`; the settings
  [`dedup()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/dedup.md)
  used; or `list(method = "precomputed")` when the input already carried
  `dedup_status`).

## Value

An S7 object with properties `results` and `metadata`. Every
construction and every edit is validated.
