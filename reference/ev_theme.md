# enrichVolcano theme and palette

Default ggplot theme + locked palette hex codes for volcano + ring
plots. Lock prevents downstream palette regressions when ggplot2
updates.

## Usage

``` r
ev_theme(base_size = 11, palette = "default")
```

## Arguments

- base_size:

  Numeric base font size.

- palette:

  Character; named palette set ("default" only for v1).

## Value

A list with `base_size`, `theme`, and `palette` elements.

## Examples

``` r
th <- ev_theme(base_size = 11)
th$palette$up
#> [1] "#D6604D"
```
