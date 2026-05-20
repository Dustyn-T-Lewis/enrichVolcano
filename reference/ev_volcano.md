# Build a volcano ggplot for one contrast

Build a volcano ggplot for one contrast

## Usage

``` r
ev_volcano(
  data,
  contrast,
  p_method = c("pi_eq2", "pi_eq1", "raw_p", "adj_p"),
  p_threshold = 0.05,
  logfc_threshold = log2(1.5),
  label_n = 10,
  label_by = NULL,
  label_genes = NULL,
  theme = ev_theme()
)
```

## Arguments

- data:

  Tidy long tibble (post `pi_score`/`adjust_p`).

- contrast:

  Single contrast name.

- p_method:

  Which p-column to use on the y-axis: `"pi_eq2"`, `"pi_eq1"`,
  `"raw_p"`, or `"adj_p"`.

- p_threshold:

  Significance cutoff (default 0.05).

- logfc_threshold:

  Effect-size cutoff (default log2(1.5)).

- label_n:

  Top N proteins to label.

- label_by:

  Column used to rank labels.

- label_genes:

  Optional character vector forcing specific labels.

- theme:

  Output of
  [`ev_theme()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/ev_theme.md).

## Value

ggplot object.
