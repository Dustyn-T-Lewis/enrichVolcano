#' Theme + palette for `volcano_ring()`
#'
#' Returns a list with a `theme` element (ggplot2 theme additions) and a
#' `palette` element giving the up / down / non-significant point colours and
#' the diverging NES ramp.
#'
#' @param base_size Numeric base font size.
#' @param base_family Character base font family; `""` defers to ggplot2.
#' @param palette One of `"default"` (red-blue diverging, the YvO 2026 lock),
#'   `"viridis"` (5-stop magma cuts), or `"okabe"` (Okabe-Ito CB-safe pair).
#'   Sets the starting up / down / non-significant colours and NES ramp; any
#'   of `up`, `down`, `ns`, `nes_colors` below override it.
#' @param up,down,ns Optional single colours overriding the palette's
#'   up-regulated, down-regulated, and non-significant point colours.
#' @param nes_colors Optional colour vector overriding the diverging NES ramp.
#'   When supplied without `nes_stops`, stops spread evenly across
#'   `nes_limits` (or `c(-3, 3)`).
#' @param nes_limits Optional length-2 numeric. When `NULL`, the colour scale
#'   in `volcano_ring()` uses its own default.
#' @param nes_stops Optional numeric vector matching the NES ramp length,
#'   overriding the palette's `nes_values`.
#' @return A list `list(base_size, base_family, palette, nes_limits)` consumed
#'   by [volcano_ring()].
#' @export
#' @examples
#' th <- volcano_ring_theme(base_size = 11)
#' th$palette$up
#'
#' # Custom point and arc colours, no list-poking needed:
#' th <- volcano_ring_theme(up = "#B2182B", down = "#2166AC", ns = "grey80")
volcano_ring_theme <- function(base_size = 11,
                               base_family = "",
                               palette = c("default", "viridis", "okabe"),
                               up = NULL,
                               down = NULL,
                               ns = NULL,
                               nes_colors = NULL,
                               nes_limits = NULL,
                               nes_stops = NULL) {
  palette <- match.arg(palette)
  for (col in list(up, down, ns)) ev_assert_colour(col, "colour override")
  palettes <- list(
    default = list(
      up = "#D6604D", down = "#4393C3", ns = "grey70",
      nes_scale = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
      nes_values = c(-3, -1.5, 0, 1.5, 3)
    ),
    viridis = list(
      up = "#FDE725", down = "#440154", ns = "grey70",
      nes_scale = c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"),
      nes_values = c(-3, -1.5, 0, 1.5, 3)
    ),
    okabe = list(
      up = "#D55E00", down = "#0072B2", ns = "grey70",
      nes_scale = c("#0072B2", "#56B4E9", "white", "#E69F00", "#D55E00"),
      nes_values = c(-3, -1.5, 0, 1.5, 3)
    )
  )
  pal <- palettes[[palette]]
  if (!is.null(up)) pal$up <- up
  if (!is.null(down)) pal$down <- down
  if (!is.null(ns)) pal$ns <- ns
  if (!is.null(nes_colors)) {
    pal$nes_scale <- nes_colors
    span <- nes_limits %||% c(-3, 3)
    pal$nes_values <- seq(span[1], span[2], length.out = length(nes_colors))
  }
  if (!is.null(nes_stops)) {
    if (length(nes_stops) != length(pal$nes_scale)) {
      ev_abort(
        c("`nes_stops` length must match the palette ramp.",
          "i" = "Expected {length(pal$nes_scale)} stops, got {length(nes_stops)}."
        ),
        class = "enrichVolcano_param_error"
      )
    }
    pal$nes_values <- nes_stops
  }
  list(
    base_size = base_size,
    base_family = base_family,
    palette = pal,
    nes_limits = nes_limits
  )
}
