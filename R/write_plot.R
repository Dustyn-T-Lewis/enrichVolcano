#' Write panels to one PDF
#'
#' Page 1 lays every panel out as a composite lettered A, B, C in list order.
#' Each later page holds one panel at full size, headed by its letter and
#' name. Fonts are embedded, so the file looks the same on any machine: macOS
#' writes through the Quartz PDF device, other systems through [grDevices::cairo_pdf()].
#'
#' @param panels A named list of ggplots, such as the output of
#'   [plot_volcano_ring()] and [plot_scatter()]. Names head the single-panel
#'   pages.
#' @param file Path of the PDF to write.
#' @param ncol,nrow Grid of the composite, passed to [patchwork::wrap_plots()].
#' @param design A patchwork layout such as `"AB\nCC"`, used instead of `ncol`
#'   and `nrow`.
#' @param caption Text under the composite. `NULL` draws none.
#' @param width,height Page size in inches. The default is US Letter landscape.
#' @return `file`, invisibly.
#' @export
#' @examples
#' da <- read_example("yvo")$da
#' ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
#'   package = "enrichVolcano"
#' )))
#' panels <- list(
#'   "Young, trained" = plot_volcano_ring(da, ex, contrast = "Training_Young"),
#'   "Old, trained" = plot_volcano_ring(da, ex, contrast = "Training_Old")
#' )
#' file <- tempfile(fileext = ".pdf")
#' write_plot(panels, file, ncol = 2)
write_plot <- function(panels, file, ncol = NULL, nrow = NULL, design = NULL, caption = NULL,
                       width = 11, height = 8.5) {
  composite <- compose_panels(panels, ncol, nrow, design, caption)
  if (capabilities("aqua")) {
    # quartz is not exported on Windows, so look it up only where it exists.
    getExportedValue("grDevices", "quartz")(type = "pdf", file = file, width = width, height = height)
  } else {
    grDevices::cairo_pdf(file, width = width, height = height, onefile = TRUE)
  }
  on.exit(grDevices::dev.off())
  print(composite)
  for (i in seq_along(panels)) {
    heading <- patchwork::plot_annotation(
      title = paste(LETTERS[i], names(panels)[i]),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    )
    print(panels[[i]] + heading)
  }
  invisible(file)
}

compose_panels <- function(panels, ncol = NULL, nrow = NULL, design = NULL, caption = NULL) {
  ok <- is.list(panels) && !inherits(panels, "ggplot") && length(panels) %in% 1:26 &&
    !is.null(names(panels)) && all(nzchar(names(panels))) &&
    all(vapply(panels, inherits, logical(1), "ggplot"))
  if (!ok) {
    ev_abort("{.arg panels} must be a named list of 1 to 26 ggplots.", class = "enrichVolcano_param_error")
  }
  patchwork::wrap_plots(panels, ncol = ncol, nrow = nrow, design = design) +
    patchwork::plot_annotation(tag_levels = "A", caption = caption)
}
