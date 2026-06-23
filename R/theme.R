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
#' @param nes_limits Optional length-2 numeric. When `NULL`, the colour scale
#'   in `volcano_ring()` uses its own default.
#' @param nes_stops Optional numeric vector matching `palette`'s `nes_scale`
#'   length, overriding the palette's `nes_values`.
#' @return A list `list(base_size, base_family, theme, palette)`.
#' @export
#' @examples
#' th <- volcano_ring_theme(base_size = 11)
#' th$palette$up
volcano_ring_theme <- function(base_size = 11,
                               base_family = "",
                               palette = c("default", "viridis", "okabe"),
                               nes_limits = NULL,
                               nes_stops = NULL) {
  palette <- match.arg(palette)
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
    theme = ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position  = "bottom"
      ),
    palette = pal,
    nes_limits = nes_limits
  )
}
