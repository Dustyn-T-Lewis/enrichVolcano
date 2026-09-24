# Convert enrichment results into an `enrichment` object

The one entry point for results from any enrichment tool. Hand it one
table per contrast and it returns a validated `enrichment` object that
every plot in the package reads.

## Usage

``` r
as_enrichment(
  x,
  enrichment_test = NULL,
  database = NULL,
  term = "term",
  score = "score",
  padj = "padj",
  p = "p",
  size = "size",
  direction = "direction",
  leading_edge = NULL,
  score_type = NULL
)
```

## Arguments

- x:

  A long table with a `contrast` column, or a named list of tables, one
  per contrast, whose names become the contrasts.

- enrichment_test:

  The method that produced `x`: one of `"fgsea"`, `"fry"`, `"mroast"`,
  `"camera"`, `"cameraPR"`, `"gseaResult"`, `"ora"`, `"custom"`. `NULL`
  recognises it from the columns.

- database:

  A label for every row, such as `"Hallmark"`, when the input has no
  `database` column.

- term, score, padj, p, size, direction, leading_edge:

  Column names, used only when `enrichment_test` is `"custom"` or
  `"ora"`. `p`, `size` and `leading_edge` are optional; a leading-edge
  column may be a list or a string separated by `;`, `/` or `|`. Pass
  `direction = NULL` to take direction from the score sign.

- score_type:

  Axis and legend label for a custom score, e.g. `"NES"`.

## Value

An S7 `enrichment` object. Every construction and every edit is
validated. `@results` has one row per term per contrast:

- `contrast`, `database`, `term`: character. `database` may be `NA`.

- `score`: numeric. NES for fgsea and clusterProfiler; signed
  \\-\log\_{10}\\(FDR) for limma rotation/competitive tests and ORA.

- `p`, `padj`: numeric in \[0, 1\]; `NA` allowed.

- `size`: numeric set size.

- `direction`: `"up"` or `"down"`, agreeing with the sign of `score`.

- `leading_edge`: list of character vectors, empty when the test has
  none.

Extra columns, such as `dedup_status`, are kept. `@metadata` holds
`enrichment_test`, `score_type` (the axis and legend label for `score`)
and `dedup` (the settings
[`dedup_terms()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/dedup_terms.md)
used, or `list(method = "precomputed")` when the input already carried
`dedup_status`).

## Which test produced the table

The test is recognised from the columns it returns:

- **fgsea**: `pathway, pval, padj, NES, size, leadingEdge`. Score is
  NES.

- **fry** (`PValue.Mixed`), **mroast** (`PropUp`), **camera** estimating
  its own inter-gene correlation (`Correlation`):
  `NGenes, Direction, PValue, FDR`, set names in the row names (or a
  `term` column, which a long table with a `contrast` column must have).
  None of these report an effect size or a leading edge, so the score is
  signed \\-\log\_{10}\\(FDR) and leading edges are empty.

- **camera** with its default fixed correlation and **cameraPR** return
  identical columns. Say which with `enrichment_test = "camera"` or
  `"cameraPR"`.

- **gseaResult** (clusterProfiler `GSEA()`, `gseGO()`, ...): read from
  the object's `result` table, or from that table exported as a data
  frame. Score is NES; leading edges come from `core_enrichment`.

- **ora** (enrichR, `enrichGO()`, g:Profiler, ...): over-representation
  has no direction of its own, so supply a `direction` column (run up-
  and down-regulated genes separately and stack them). Score is signed
  \\-\log\_{10}\\(padj).

- **custom**: any table, mapped with the column arguments. Direction
  comes from the `direction` column when present, otherwise from the
  score sign.

fry and mroast are self-contained tests (is this set changed at all?);
camera, cameraPR and fgsea are competitive (is it changed more than the
genes outside it?). The test is recorded so figures can say which.

## References

Wu T, Hu E, Xu S, et al. (2021). clusterProfiler 4.0: A universal
enrichment tool for interpreting omics data. The Innovation 2(3):100141.
[doi:10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141)

## Examples

``` r
path <- system.file("extdata", "examples", "yvo_fgsea.csv.gz",
  package = "enrichVolcano"
)
ex <- as_enrichment(read.csv(path))
#> Reading "fgsea" results.
ex
#> <enrichment> fgsea, score: NES
#> 4 contrasts: Aging, Training_Young, Training_Old, and Interaction
#> 5570 rows: GO Slim (108), GO:BP (3946), Hallmark (144), KEGG (84), and Reactome
#> (1288)
#> Dedup: precomputed
```
