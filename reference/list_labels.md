# List the term labels a figure will draw, to review or edit

Picks terms the way
[`plot_volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_volcano_ring.md)
does, for each contrast, and lists each term once with both shipped
styles from
[`clean_label()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/clean_label.md).
Edit the `label` column, in R or in the CSV that `write_labels()`
writes, and pass the table back through `labels` of
[`plot_volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_volcano_ring.md)
or
[`plot_scatter()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/plot_scatter.md).
A line break is written as the two characters `\n`.

## Usage

``` r
list_labels(
  enrichment,
  contrast = NULL,
  databases = NULL,
  collapse = TRUE,
  term_threshold = 0.05,
  n_terms = 12
)

write_labels(enrichment, file, ...)
```

## Arguments

- enrichment:

  An enrichment object from
  [`as_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/as_enrichment.md)
  or
  [`run_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/run_enrichment.md).

- contrast:

  Contrasts to list; `NULL` (default) lists every contrast.

- databases:

  Collections to draw from, matched against the `database` column, e.g.
  `c("Hallmark", "GO Slim")`. `NULL` (default) draws from all. Skipped,
  with a note, when the object has no database labels.

- collapse:

  Hide terms flagged `"redundant"` by
  [`dedup_terms()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/dedup_terms.md).

- term_threshold:

  Terms need `padj` below this to be drawn.

- n_terms:

  Most unique terms per contrast, as on the ring. `Inf` lists every
  significant term.

- file:

  Path of the CSV to write.

- ...:

  Arguments passed on to `list_labels()`.

## Value

A data frame with one row per term: `database`, `term`, `contrasts` (the
contrasts that draw it, joined by `;`), `clean`, `short` and `label`, a
copy of `short` to edit.

`write_labels()`: `file`, invisibly.

## Examples

``` r
ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
  package = "enrichVolcano"
)))
#> Reading "fgsea" results.
head(list_labels(ex))
#>   database                                                  term
#> 1 Reactome                                  REACTOME_TRANSLATION
#> 2 Reactome               REACTOME_RESPIRATORY_ELECTRON_TRANSPORT
#> 3     KEGG KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_I
#> 4 Reactome          REACTOME_RIBOSOME_ASSOCIATED_QUALITY_CONTROL
#> 5    GO:BP         GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT
#> 6  GO Slim                        GOSLIM_CYTOPLASMIC_TRANSLATION
#>                    contrasts                                    clean
#> 1         Aging;Training_Old                              Translation
#> 2 Training_Young;Interaction           Respiratory Electron Transport
#> 3             Training_Young Reference Electron Transfer In Complex I
#> 4                      Aging      Ribosome Associated Quality Control
#> 5 Training_Young;Interaction ATP Synthesis Coupled Electron Transport
#> 6         Aging;Training_Old                  Cytoplasmic Translation
#>                                      short
#> 1                              Translation
#> 2           Respiratory Electron Transport
#> 3 Reference Electron Transfer In Complex I
#> 4      Ribosome Associated Quality Control
#> 5 ATP Synthesis Coupled Electron Transport
#> 6                  Cytoplasmic Translation
#>                                      label
#> 1                              Translation
#> 2           Respiratory Electron Transport
#> 3 Reference Electron Transfer In Complex I
#> 4      Ribosome Associated Quality Control
#> 5 ATP Synthesis Coupled Electron Transport
#> 6                  Cytoplasmic Translation
```
