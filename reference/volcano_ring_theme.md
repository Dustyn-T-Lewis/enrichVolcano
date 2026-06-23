# Theme + palette for `volcano_ring()`

Returns a list with a `theme` element (ggplot2 theme additions) and a
`palette` element giving the up / down / non-significant point colours
and the diverging NES ramp.

## Usage

``` r
volcano_ring_theme(
  base_size = 11,
  base_family = "",
  palette = c("default", "viridis", "okabe"),
  nes_limits = NULL,
  nes_stops = NULL
)
```

## Arguments

- base_size:

  Numeric base font size.

- base_family:

  Character base font family; `""` defers to ggplot2.

- palette:

  One of `"default"` (red-blue diverging, the YvO 2026 lock),
  `"viridis"` (5-stop magma cuts), or `"okabe"` (Okabe-Ito CB-safe
  pair).

- nes_limits:

  Optional length-2 numeric. When `NULL`, the colour scale in
  [`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
  uses its own default.

- nes_stops:

  Optional numeric vector matching `palette`'s `nes_scale` length,
  overriding the palette's `nes_values`.

## Value

A list `list(base_size, base_family, theme, palette)`.

## Examples

``` r
th <- volcano_ring_theme(base_size = 11)
th$palette$up
#> [1] "#D6604D"
```
