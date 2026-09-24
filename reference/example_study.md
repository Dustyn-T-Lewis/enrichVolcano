# Example studies shipped with the package

Seven proteomics studies from one lab, prepared as
[`read_study()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/read_study.md)
inputs. The three limpa studies (`bfr_limpa`, `mouse_pas`, `yvo`)
include the matrix, samples, contrasts and precision weights, so every
test in
[`run_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/run_enrichment.md)
runs on them; the rest ship DA results only, for fgsea and the figures.
Call with no name to list them with their species and designs.

## Usage

``` r
example_study(name = NULL)
```

## Arguments

- name:

  A study name, or `NULL` to list them.

## Value

The study, as
[`read_study()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/read_study.md)
returns it, with its index row in `info`; or the index as a data frame.

## Examples

``` r
example_study()
#>           name           species subject_effect
#> 1    bfr_limpa      Homo sapiens          fixed
#> 2    mouse_pas      Mus musculus           <NA>
#> 3          yvo      Homo sapiens           <NA>
#> 4 bfr_proteoda      Homo sapiens           <NA>
#> 5          cvh      Homo sapiens           <NA>
#> 6        hrvlr      Homo sapiens           <NA>
#> 7         mito Rattus norvegicus           <NA>
#>                                                                                        description
#> 1 Human muscle, one leg blood-flow-restricted and one heavy-load, before and after training; limpa
#> 2                                               Mouse muscle mitochondria, treatment by sex; limpa
#> 3                                    Human muscle, young and old, before and after training; limpa
#> 4                                         The blood-flow-restriction cohort analysed with proteoDA
#> 5                            Human muscle, cancer recovery arms against healthy controls; proteoDA
#> 6                          Human muscle, high and low responders across three timepoints; proteoDA
#> 7                              Rat cardiomyoblasts, stressor by mitochondrial transplant; proteoDA
```
