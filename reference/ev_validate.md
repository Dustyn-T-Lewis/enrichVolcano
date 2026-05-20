# Validate and canonicalize input from any supported source

Detects input format from column patterns, applies the appropriate
adapter, and returns a tidy long tibble with canonical columns (gene,
contrast, logFC, P.Value, adj.P.Val).

## Usage

``` r
ev_validate(data, x = NULL, y = NULL, lab = NULL)
```

## Arguments

- data:

  Input data frame from any supported source.

- x:

  Column name for log fold change (defaults to auto-detect).

- y:

  Column name for p-value (defaults to auto-detect).

- lab:

  Column name for gene labels (defaults to auto-detect).

## Value

Tidy long tibble with attribute `ev_source` indicating detected format.
