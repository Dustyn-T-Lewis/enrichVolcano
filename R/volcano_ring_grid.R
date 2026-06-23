#' Compose a grid of `volcano_ring()` plots, one per contrast
#'
#' @param volc_dfs Named list of tidy DA tibbles, one per contrast. A single
#'   data.frame carrying a `contrast` column is also accepted and is split.
#' @param enrich_dfs Named list of tidy enrichment tibbles, one per contrast.
#'   Same split convenience as `volc_dfs`.
#' @param contrasts Character vector of contrast names to include and the
#'   order in which to draw them. Defaults to `names(volc_dfs)`.
#' @param nrow,ncol Outer layout dims; forwarded to `patchwork::wrap_plots()`.
#' @param tag_levels Panel-tag scheme; forwarded to `patchwork::plot_annotation()`.
#' @param guides Patchwork `guides` argument; default `"collect"` collects the
#'   shared NES legend.
#' @param ... Forwarded to each `volcano_ring()` call (e.g. `gene_col`,
#'   `padj_col`, `magnitude`, `theme`).
#' @return An S3 object `c("volcano_ring_grid", "list")` with elements
#'   `$plot` (patchwork) and `$data` (list of `list(volc, enrich)` pairs).
#' @export
volcano_ring_grid <- function(volc_dfs, enrich_dfs,
                              contrasts = NULL,
                              nrow = NULL,
                              ncol = NULL,
                              tag_levels = "A",
                              guides = "collect",
                              ...) {
  volc_list <- split_by_contrast(volc_dfs, "volc_dfs")
  enrich_list <- split_by_contrast(enrich_dfs, "enrich_dfs")

  contrasts <- contrasts %||% names(volc_list)
  if (is.null(contrasts) || !length(contrasts)) {
    ev_abort(
      "Cannot infer contrasts; supply `contrasts` or name the list elements.",
      class = "enrichVolcano_input_error"
    )
  }
  missing_v <- setdiff(contrasts, names(volc_list))
  missing_e <- setdiff(contrasts, names(enrich_list))
  if (length(missing_v) || length(missing_e)) {
    ev_abort(
      c("Some contrasts are missing input frames.",
        "i" = "Missing in volc_dfs: {.val {missing_v}}",
        "i" = "Missing in enrich_dfs: {.val {missing_e}}"
      ),
      class = "enrichVolcano_input_error"
    )
  }

  panels <- lapply(contrasts, function(cn) {
    volcano_ring(volc_list[[cn]], enrich_list[[cn]], title = cn, ...)
  })
  names(panels) <- contrasts

  plot <- patchwork::wrap_plots(panels, nrow = nrow, ncol = ncol, guides = guides) +
    patchwork::plot_annotation(tag_levels = tag_levels)

  out <- list(
    plot = plot,
    data = stats::setNames(
      lapply(contrasts, function(cn) {
        list(
          volc = volc_list[[cn]],
          enrich = enrich_list[[cn]]
        )
      }),
      contrasts
    )
  )
  class(out) <- c("volcano_ring_grid", "list")
  out
}

#' @export
print.volcano_ring_grid <- function(x, ...) {
  print(x$plot, ...)
  invisible(x)
}

split_by_contrast <- function(x, arg_name) {
  if (is.list(x) && !is.data.frame(x)) {
    return(x)
  }
  if (is.data.frame(x)) {
    if (!"contrast" %in% names(x)) {
      ev_abort(
        c("{.arg {arg_name}} is a single data.frame but has no `contrast` column.",
          "i" = "Either pass a named list of frames, or add a `contrast` column."
        ),
        class = "enrichVolcano_input_error"
      )
    }
    return(split(x, x$contrast))
  }
  ev_abort(
    "{.arg {arg_name}} must be a named list of data.frames or a single data.frame.",
    class = "enrichVolcano_input_error"
  )
}
