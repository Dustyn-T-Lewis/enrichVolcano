# Compose volcano and ring panels into a patchwork

Compose volcano and ring panels into a patchwork

## Usage

``` r
ev_compose(
  volcano_plots,
  ring_plots,
  nrow = NULL,
  ncol = NULL,
  guides = "collect",
  data = NULL
)
```

## Arguments

- volcano_plots:

  Named list of ggplots (one per contrast).

- ring_plots:

  Named list of ggplots (one per contrast).

- nrow, ncol:

  Outer layout dims (passed to
  [`patchwork::wrap_plots`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)).

- guides:

  Patchwork guides argument; default `"collect"`.

- data:

  Optional list of intermediate tibbles to attach as `ev_data`.

## Value

Object of class `c("enrichVolcano", "patchwork", "gg", "ggplot")`.
