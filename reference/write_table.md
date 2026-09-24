# Write enrichment results to a CSV file

Writes one row per term per contrast, with an `enrichment_test` column
and leading edges joined by `;`.
[`as_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_enrichment.md)
reads the file back with `enrichment_test = "custom"` and
`leading_edge = "leading_edge"`.

## Usage

``` r
write_table(x, file)
```

## Arguments

- x:

  An enrichment object, or the list of them that
  [`run_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/run_enrichment.md)
  returns, which is written as one long table.

- file:

  Path of the CSV to write.

## Value

`file`, invisibly.

## Examples

``` r
ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
  package = "enrichVolcano"
)))
#> Reading "fgsea" results.
write_table(ex, tempfile(fileext = ".csv"))
```
