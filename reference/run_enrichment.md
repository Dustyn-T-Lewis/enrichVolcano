# Run gene-set enrichment on a study

Tests each contrast of a study against gene-set collections and returns
one enrichment object per test, the same object
[`as_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_enrichment.md)
builds.

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

  Test sets with at least `min_size` and at most `max_size` genes
  present in the data.

## Value

A named list with one enrichment object per test. Its metadata records
the ranking statistic and the gene-set versions.

## fgsea

Proteins are ranked by the study's `rank` column (see
[`as_da()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_da.md)).
Where several proteins share a gene symbol, the most abundant one
represents the gene (the table's average abundance, else the matrix row
mean); without any abundance, the first accession alphabetically. The
choice never looks at the results. Each collection is tested and
corrected separately. fgsea permutes, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible p-values. Each contrast and collection draws its own seed
from it, so adding or removing a collection leaves the others' results
unchanged.

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

## References

Subramanian A, Tamayo P, Mootha VK, et al. (2005). Gene set enrichment
analysis: a knowledge-based approach for interpreting genome-wide
expression profiles. PNAS 102(43):15545-15550.
[doi:10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)

Korotkevich G, Sukhov V, Sergushichev A (2019). Fast gene set enrichment
analysis. bioRxiv. [doi:10.1101/060012](https://doi.org/10.1101/060012)

Wu D, Smyth GK (2012). Camera: a competitive gene set test accounting
for inter-gene correlation. Nucleic Acids Research 40(17):e133.
[doi:10.1093/nar/gks461](https://doi.org/10.1093/nar/gks461)

Wu D, Lim E, Vaillant F, et al. (2010). ROAST: rotation gene set tests
for complex microarray experiments. Bioinformatics 26(17):2176-2182.
[doi:10.1093/bioinformatics/btq401](https://doi.org/10.1093/bioinformatics/btq401)

Smyth GK, Michaud J, Scott HS (2005). Use of within-array replicate
spots for assessing differential expression in microarray experiments.
Bioinformatics 21(9):2067-2075.
[doi:10.1093/bioinformatics/bti270](https://doi.org/10.1093/bioinformatics/bti270)

Ritchie ME, Phipson B, Wu D, et al. (2015). limma powers differential
expression analyses for RNA-sequencing and microarray studies. Nucleic
Acids Research 43(7):e47.
[doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)
