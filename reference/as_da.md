# Standardise differential-abundance results

Reads a DA results table from any common proteomics tool and returns it
with one set of column names, so the plots and
[`run_enrichment()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/run_enrichment.md)
can read it.

## Usage

``` r
as_da(
  x,
  contrast = NULL,
  species = "Homo sapiens",
  protein = NULL,
  gene = NULL,
  logfc = NULL,
  t = NULL,
  p = NULL,
  padj = NULL,
  abundance = NULL
)
```

## Arguments

- x:

  A results table with a contrast column (`contrast` or MSstats'
  `Label`), a single-contrast table with `contrast` set, or a named list
  of tables whose names are the contrasts.

- contrast:

  Label for a table that has no contrast column.

- species:

  For UniProt accessions, the species whose annotation package supplies
  current gene symbols (`"Homo sapiens"`, `"Mus musculus"`,
  `"Rattus norvegicus"`). `NULL` keeps the table's own symbols. Tables
  keyed by symbols rather than accessions keep theirs either way.

- protein, gene, logfc, t, p, padj, abundance:

  Column names, when detection does not find them. `logfc` names the
  fold-change column.

## Value

A data frame of class `enrichVolcano_da` with columns `protein`, `gene`,
`contrast`, `logFC`, `t`, `p`, `padj`, `abundance`, `rank` and
`rank_stat`, followed by any input columns it did not use (for example a
pi-value column for `volcano_ring(volc_sig_col = )`).

## Recognised tools

limma and limpa `topTable()` (IDs in the row names, or an unnamed first
column once written to CSV), proteoDA results and per-contrast CSVs,
MSstats `groupComparison()`, proDA `test_diff()`, msqrob2
`topFeatures()`, prolfQua contrast tables and ProtRank output. Any other
table works through the column arguments. F-test tables (no fold change)
and wide tables with one column per contrast are refused; reshape them
to one row per protein and contrast first.

## Ranking statistic

`rank` is the moderated t when the table has one, else
`sign(logFC) * -log10(p)`, else `sign(logFC) * -log10(padj)`;
`rank_stat` says which.

## Gene symbols

Gene sets list gene symbols, so each UniProt accession is mapped to its
current symbol through the species' annotation package (isoform suffixes
such as `-2` are dropped first). A symbol from the search engine's FASTA
can be out of date: `O00483` is `COXFA4`, formerly `NDUFA4`.

## Examples

``` r
tbl <- data.frame(
  logFC = c(1.2, -0.8), t = c(4.1, -3.2), P.Value = c(1e-4, 2e-3),
  adj.P.Val = c(1e-3, 0.01), row.names = c("P31040", "Q9UBK2")
)
as_da(tbl, contrast = "Aging")
#> 
#> 2 of 2 accessions mapped to Homo sapiens symbols.
#>   protein     gene contrast logFC    t     p  padj abundance rank rank_stat
#> 1  P31040     SDHA    Aging   1.2  4.1 1e-04 0.001        NA  4.1         t
#> 2  Q9UBK2 PPARGC1A    Aging  -0.8 -3.2 2e-03 0.010        NA -3.2         t
```
