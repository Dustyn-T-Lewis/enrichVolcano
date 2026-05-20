# Adjust p-values

Adjust p-values

## Usage

``` r
adjust_p(
  data,
  method = c("BH", "bonferroni", "qvalue", "IHW"),
  covariate = NULL,
  alpha = 0.05,
  p_col = "P.Value",
  by_contrast = TRUE
)
```

## Arguments

- data:

  Tidy long tibble.

- method:

  One of "BH", "bonferroni", "qvalue", "IHW".

- covariate:

  Column name (for IHW).

- alpha:

  Significance threshold (for IHW).

- p_col:

  Raw p-value column name.

- by_contrast:

  Logical; adjust within each contrast (default TRUE).

## Value

Tibble with `adj.P.Val` column (added or overwritten).

## Details

- **BH** (Benjamini-Hochberg, default): field-standard for proteomics.

- **Bonferroni**: strict family-wise; use for confirmation.

- **q-value** (Storey 2003): estimates pi_0, more power than BH;
  requires `qvalue` package.

- **IHW** (Ignatiadis 2016): BH weighted by covariate; requires `IHW`
  package and a numeric `covariate` column.
