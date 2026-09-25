#' Write panels to one PDF
#'
#' Page 1 lays every panel out as a composite lettered A, B, C in list order.
#' Each later page holds one panel at full size, headed by its letter and
#' name. [grDevices::cairo_pdf()] embeds the fonts, so the file looks the same
#' on any machine. On macOS, R's cairo needs XQuartz (<https://www.xquartz.org>).
#'
#' @param panels A named list of ggplots, such as the output of
#'   [plot_volcano_ring()] and [plot_scatter()]. Names head the single-panel
#'   pages.
#' @param file Path of the PDF to write.
#' @param ncol,nrow Grid of the composite, passed to [patchwork::wrap_plots()].
#' @param design A patchwork layout such as `"AB\nCC"`, used instead of `ncol`
#'   and `nrow`.
#' @param caption Text under the composite. `NULL` draws none.
#' @param width,height Page size in inches. `NULL` gives each composite cell
#'   7 by 7.5 inches, the size the panels are drawn for, and widens the page to
#'   at least landscape proportions. Every page has this size.
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
#' try(write_plot(panels, file, ncol = 2))
write_plot <- function(panels, file, ncol = NULL, nrow = NULL, design = NULL, caption = NULL,
                       width = NULL, height = NULL) {
  composite <- compose_panels(panels, ncol, nrow, design, caption)
  size <- page_size(length(panels), ncol, nrow, design)
  width <- width %||% size[1]
  height <- height %||% size[2]
  before <- grDevices::dev.cur()
  # A missing cairo library only warns and opens no device, so compare devices.
  suppressWarnings(grDevices::cairo_pdf(file, width = width, height = height, onefile = TRUE))
  if (identical(grDevices::dev.cur(), before)) {
    ev_abort(
      c(
        "cairo could not start, so fonts cannot be embedded.",
        i = "On macOS, install XQuartz: {.url https://www.xquartz.org}."
      ),
      class = "enrichVolcano_device_error"
    )
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

# Page size in inches: 7 x 7.5 per composite cell, widened to at least the
# 11:8.5 landscape ratio. Mirrors patchwork's default grid.
page_size <- function(n, ncol = NULL, nrow = NULL, design = NULL) {
  if (!is.null(design)) {
    rows <- Filter(nzchar, trimws(strsplit(design, "\n")[[1]]))
    nrow <- length(rows)
    ncol <- max(nchar(rows))
  }
  if (is.null(ncol) && is.null(nrow)) ncol <- ceiling(sqrt(n))
  if (is.null(nrow)) nrow <- ceiling(n / ncol)
  if (is.null(ncol)) ncol <- ceiling(n / nrow)
  height <- 7.5 * nrow
  c(max(7 * ncol, height * 11 / 8.5), height)
}
