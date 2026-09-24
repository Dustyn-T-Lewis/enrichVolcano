#' Every significant term of one contrast around its volcano
#'
#' Draws all significant terms of a contrast as a ring, without names, to show
#' whether the contrast leans up or down. Every arc gets the same angle, so the
#' up and down halves grow with their term counts. Fill and arc height are the
#' score over the strongest drawn score: the strongest term is full colour and
#' full height, and every other term is a fraction of it.
#'
#' @inheritParams plot_volcano_ring
#' @param subtitle Subtitle. `NULL` reports the up and down term counts.
#' @param ... Other [plot_volcano_ring()] arguments that set the layout, such
#'   as `ring_radius`, `arc_height_range`, `p_threshold` or `point_size`.
#' @return A ggplot.
#' @export
#' @examples
#' da <- read_example("yvo")$da
#' ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
#'   package = "enrichVolcano"
#' )))
#' plot_bias_ring(da, ex, contrast = "Aging", title = "Aging")
plot_bias_ring <- function(da, enrichment, contrast = NULL, databases = NULL, collapse = TRUE,
                           term_threshold = 0.05, title = NULL, subtitle = NULL,
                           theme = plot_theme(), ...) {
  defaults <- formals(plot_volcano_ring)[-(1:2)]
  args <- lapply(defaults, eval, envir = environment(plot_volcano_ring))
  extra <- list(...)
  unknown <- setdiff(names(extra), names(defaults))
  if (length(unknown) > 0 || any(!nzchar(names(extra)))) {
    ev_abort("Unknown argument{?s} in {.arg ...}: {.val {unknown}}.", class = "enrichVolcano_param_error")
  }
  given <- list(
    da = da, enrichment = enrichment, contrast = contrast, databases = databases,
    collapse = collapse, term_threshold = term_threshold, title = title,
    subtitle = subtitle, theme = theme, n_terms = Inf
  )
  args[names(extra)] <- extra
  args[names(given)] <- given
  draw_ring(args, term_labels = FALSE, normalise = TRUE)
}
