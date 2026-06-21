#' enrichVolcano theme and palette
#'
#' Default ggplot theme + locked palette hex codes for volcano + ring plots.
#' Lock prevents downstream palette regressions when ggplot2 updates.
#'
#' @param base_size Numeric base font size.
#' @param palette Character; named palette set ("default" only for v1).
#' @return A list with `base_size`, `theme`, and `palette` elements.
#' @export
#' @examples
#' th <- ev_theme(base_size = 11)
#' th$palette$up
ev_theme <- function(base_size = 11, palette = "default") {
  palettes <- list(
    default = list(
      up = "#D6604D",
      down = "#4393C3",
      ns = "grey70",
      # Single source of truth for the NES diverging ramp.
      # `ring_plot()` uses indices 1, 3, 5 in scale_fill_gradient2();
      # `ev_volcano_ring()` uses the full 5 stops in scale_fill_gradientn().
      nes_scale  = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
      nes_values = c(-3, -1.5, 0, 1.5, 3)
    )
  )
  if (!palette %in% names(palettes)) {
    ev_abort(
      "Unknown palette {.val {palette}}. Available: {.val {names(palettes)}}.",
      class = "ev_unknown_palette"
    )
  }
  list(
    base_size = base_size,
    theme = ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom"
      ),
    palette = palettes[[palette]]
  )
}
