# Run gene-set enrichment on a study

Tests each contrast of a study against gene-set collections and returns
one
[enrichment](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrichment.md)
object per test, ready for the plots.

## Usage

``` r
run_enrichment(
  study,
  gene_sets,
  tests = c("fgsea", "camera", "fry"),
  subject_effect = c("block", "fixed"),
  covariates = NULL,
  inter_gene_cor = 0.01,
  min_size = 15,
  max_size = 500
)
```

## Arguments

- study:

  Output of
  [`read_study()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/read_study.md),
  or of
  [`as_da()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_da.md)
  for fgsea alone.

- gene_sets:

  Output of
  [`load_gene_sets()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/load_gene_sets.md),
  a named list of collections of gene-symbol vectors, or one named list
  of sets.

- tests:

  Any of `"fgsea"`, `"camera"`, `"fry"`. A study without a matrix runs
  fgsea only.

- subject_effect:

  `"block"` or `"fixed"`: how a `subject` column enters the camera and
  fry model.

- covariates:

  Sample columns added to the camera and fry design.

- inter_gene_cor:

  camera's inter-gene correlation: a number, or `NA` to estimate it
  (unblocked designs only).

- min_size, max_size:

  Sets need this many genes present in the data.

## Value

A named list with one
[enrichment](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/enrichment.md)
object per test. Its metadata records the ranking statistic and the
gene-set versions.

## fgsea

Proteins are ranked by the study's `rank` column (see
[`as_da()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_da.md)).
Where several proteins share a gene symbol, the most abundant one
represents the gene (the table's average abundance, else the matrix row
mean); without any abundance, the first accession alphabetically. The
choice never looks at the results. Each collection is tested and
corrected separately. fgsea permutes, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible p-values.

## camera and fry

Both refit limma on the study's matrix, which must have no missing
values (use an imputed or limpa-quantified matrix), reduced to one
protein per gene by the same rule. The design is `~ 0 + group`, plus any
`covariates`, plus `subject` when `subject_effect = "fixed"`; contrasts
come from the `contrasts` sheet. With `subject_effect = "block"` and a
`subject` column, samples from one subject are treated as correlated
([`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html)).
fry uses that blocking directly. camera cannot block, so for a blocked
design the package runs
[`limma::cameraPR()`](https://rdrr.io/pkg/limma/man/camera.html) on the
moderated t of the blocked fit instead, and the result records
`cameraPR` as its test. Precision weights are used when the study has
them. camera assumes genes within a set correlate at `inter_gene_cor`
(limma's default, 0.01); `NA` estimates it from the data, which can
change results substantially and needs an unblocked design.
