#' Compose volcano and ring panels into a patchwork
#'
#' @param volcano_plots Named list of ggplots (one per contrast).
#' @param ring_plots Named list of ggplots (one per contrast).
#' @param nrow,ncol Outer layout dims (passed to `patchwork::wrap_plots`).
#' @param guides Patchwork guides argument; default `"collect"`.
#' @param data Optional list of intermediate tibbles to attach as `ev_data`.
#' @return Object of class `c("enrichVolcano", "patchwork", "gg", "ggplot")`.
#' @export
ev_compose <- function(volcano_plots, ring_plots,
                       nrow = NULL, ncol = NULL,
                       guides = "collect",
                       data = NULL) {
  if (!identical(names(volcano_plots), names(ring_plots))) {
    ev_abort("volcano_plots and ring_plots must have identical names.",
             class = "ev_compose_name_mismatch")
  }
  per_contrast <- purrr::map2(volcano_plots, ring_plots, function(v, r) {
    patchwork::wrap_plots(v, r, ncol = 1)
  })
  combined <- patchwork::wrap_plots(per_contrast,
                                    nrow = nrow, ncol = ncol,
                                    guides = guides)
  class(combined) <- c("enrichVolcano", class(combined))
  if (!is.null(data)) attr(combined, "ev_data") <- data
  combined
}

#' Print an enrichVolcano composite
#'
#' Prints the reconstructible call and a one-line data summary before
#' rendering the underlying patchwork.
#'
#' @param x An `enrichVolcano` object.
#' @param ... Passed to patchwork print method.
#' @export
print.enrichVolcano <- function(x, ...) {
  call_attr <- attr(x, "ev_call")
  data_attr <- attr(x, "ev_data")
  cli::cli_h2("enrichVolcano")
  if (!is.null(call_attr)) {
    cli::cli_text("Call: {.code {paste(deparse(call_attr), collapse = ' ')}}")
  }
  if (!is.null(data_attr$validated_input)) {
    cli::cli_text("Contrasts: {length(unique(data_attr$validated_input$contrast))}")
    cli::cli_text("Proteins: {nrow(data_attr$validated_input)}")
  }
  if (!is.null(data_attr$dedup_result)) {
    n_sig <- sum(data_attr$dedup_result$dedup_kept, na.rm = TRUE)
    cli::cli_text("Significant pathways (post-dedup): {n_sig}")
  }
  NextMethod()
}
