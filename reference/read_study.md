# Read a study from a workbook or a folder of CSV files

A study is up to five tables, linked by name:

- `da_results` (required): differential-abundance results in any format
  [`as_da()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_da.md)
  reads, with a contrast column.

- `matrix`: one row per protein; the first column holds the protein IDs
  and every other column is a sample.

- `samples`: `sample` and `group`, optionally `subject` for repeated
  measures, plus any covariate columns.

- `contrasts`: `name` and `expression`, arithmetic on group names such
  as `Post - Pre` or `(B_Post - B_Pre) - (A_Post - A_Pre)`.

- `weights` (optional): precision weights shaped like `matrix`.

## Usage

``` r
read_study(path, species = "Homo sapiens")
```

## Arguments

- path:

  An `.xlsx` workbook with sheets named as above, or a folder of
  `<sheet>.csv` or `<sheet>.csv.gz` files.

- species:

  Passed to
  [`as_da()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_da.md)
  for the gene-symbol lookup.

## Value

A list of class `enrichVolcano_study` with `da`, `matrix`, `samples`,
`contrasts` and `weights` (absent parts are `NULL`).

## Details

`da_results` alone is enough for fgsea and the figures; fry and camera
also need `matrix`, `samples` and `contrasts`. Every link is checked by
name: matrix columns against `samples$sample`, contrast expressions
against `samples$group`, `da_results` contrasts against
`contrasts$name`, and `weights` against `matrix`. Group names must be
syntactic (letters, digits, `.` and `_`), because `-` inside a name
would read as subtraction.
