#' Fonts and colours for the plots
#'
#' Sets the font and the colours that [plot_volcano_ring()] and
#' [plot_scatter()] use: the up, down and non-significant point colours, the
#' colour of terms significant in both scatter contrasts, and the diverging
#' ramp that fills arcs by score.
#'
#' @param base_size Base font size.
#' @param base_family Base font family. `""` uses the ggplot2 default.
#' @param palette `"default"` (red and blue), `"viridis"` (purple, teal and yellow) or
#'   `"okabe"` (Okabe-Ito, safe for colour-blind readers). `up`, `down`, `ns`,
#'   `both` and `score_colours` override single parts of it.
#' @param up,down,ns Colours of up-regulated, down-regulated and
#'   non-significant points. On the scatter, `up` and `down` also mark terms
#'   significant only in `x` or only in `y`, and tint the same-direction and
#'   opposite-direction quadrants.
#' @param both Scatter colour of terms significant in both contrasts.
#' @param score_colours Colours of the score ramp, low to high. Without
#'   `score_stops` they spread evenly across `score_limits`, or `c(-3, 3)`.
#' @param score_limits Length-2 limits of the score scale. `NULL` lets each
#'   plot choose; see `score_limits` in [plot_volcano_ring()].
#' @param score_stops Score values at which each ramp colour sits, one per
#'   colour.
#' @return A list with `base_size`, `base_family`, `palette` and
#'   `score_limits`, passed to the plots as `theme`.
#' @export
#' @examples
#' th <- plot_theme(base_size = 11)
#' th$palette$up
#'
#' th <- plot_theme(up = "#B2182B", down = "#2166AC", ns = "grey80")
plot_theme <- function(base_size = 11,
                       base_family = "",
                       palette = c("default", "viridis", "okabe"),
                       up = NULL,
                       down = NULL,
                       ns = NULL,
                       both = NULL,
                       score_colours = NULL,
                       score_limits = NULL,
                       score_stops = NULL) {
  palette <- match.arg(palette)
  for (col in list(up, down, ns, both)) ev_assert_colour(col, "colour override")
  palettes <- list(
    default = list(
      up = "#D6604D", down = "#4393C3", ns = "grey70", both = "#2E7D32",
      nes_scale = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
      nes_values = c(-3, -1.5, 0, 1.5, 3)
    ),
    viridis = list(
      up = "#21908C", down = "#440154", ns = "grey70", both = "#5DC863",
      nes_scale = c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"),
      nes_values = c(-3, -1.5, 0, 1.5, 3)
    ),
    okabe = list(
      up = "#D55E00", down = "#0072B2", ns = "grey70", both = "#009E73",
      nes_scale = c("#0072B2", "#56B4E9", "white", "#E69F00", "#D55E00"),
      nes_values = c(-3, -1.5, 0, 1.5, 3)
    )
  )
  pal <- palettes[[palette]]
  if (!is.null(up)) pal$up <- up
  if (!is.null(down)) pal$down <- down
  if (!is.null(ns)) pal$ns <- ns
  if (!is.null(both)) pal$both <- both
  if (!is.null(score_colours)) {
    pal$nes_scale <- score_colours
    span <- score_limits %||% c(-3, 3)
    pal$nes_values <- seq(span[1], span[2], length.out = length(score_colours))
  }
  if (!is.null(score_stops)) {
    if (length(score_stops) != length(pal$nes_scale)) {
      ev_abort(
        c("`score_stops` length must match the palette ramp.",
          "i" = "Expected {length(pal$nes_scale)} stops, got {length(score_stops)}."
        ),
        class = "enrichVolcano_param_error"
      )
    }
    pal$nes_values <- score_stops
  }
  list(
    base_size = base_size,
    base_family = base_family,
    palette = pal,
    score_limits = score_limits
  )
}
