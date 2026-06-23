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
  up = NULL,
  down = NULL,
  ns = NULL,
  nes_colors = NULL,
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
  pair). Sets the starting up / down / non-significant colours and NES
  ramp; any of `up`, `down`, `ns`, `nes_colors` below override it.

- up, down, ns:

  Optional single colours overriding the palette's up-regulated,
  down-regulated, and non-significant point colours.

- nes_colors:

  Optional colour vector overriding the diverging NES ramp. When
  supplied without `nes_stops`, stops spread evenly across `nes_limits`
  (or `c(-3, 3)`).

- nes_limits:

  Optional length-2 numeric. When `NULL`, the colour scale in
  [`volcano_ring()`](https://Dustyn-T-Lewis.github.io/enrichVolcano/reference/volcano_ring.md)
  uses its own default.

- nes_stops:

  Optional numeric vector matching the NES ramp length, overriding the
  palette's `nes_values`.

## Value

A list `list(base_size, base_family, theme, palette)`.

## Examples

``` r
th <- volcano_ring_theme(base_size = 11)
th$palette$up
#> [1] "#D6604D"

# Custom point and arc colours, no list-poking needed:
th <- volcano_ring_theme(up = "#B2182B", down = "#2166AC", ns = "grey80")
```
