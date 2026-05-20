# Compute pi-score (Xiao 2014)

Compute pi-score (Xiao 2014)

## Usage

``` r
pi_score(
  data,
  variant = c("eq2", "eq1"),
  p_col = "P.Value",
  logfc_col = "logFC",
  out_col = NULL
)
```

## Arguments

- data:

  Tidy long tibble with `logFC` and `P.Value` columns.

- variant:

  `"eq2"` (default) or `"eq1"`.

- p_col:

  P-value column name.

- logfc_col:

  Log fold-change column name.

- out_col:

  Output column name (defaults to `pi_eq2` / `pi_eq1`).

## Value

Input tibble with one new column.

## Details

Two variants are supported, both from Xiao et al. 2014 *Bioinformatics*
30(6):801:

- **Eq.2 (default):** \\\pi = P^{\|logFC\|}\\. Range \[0,1\], smaller =
  more significant; convention threshold 0.05.

- **Eq.1 (alternative):** \\\pi = \|log_2 FC\| \times -log\_{10} P\\.
  Unbounded, larger = more significant.

## References

Xiao Y et al. (2014) *Bioinformatics* 30(6):801–807.
