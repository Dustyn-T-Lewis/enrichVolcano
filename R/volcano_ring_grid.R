#' Compose a grid of `volcano_ring()` plots, one per contrast
#'
#' @param volc_dfs Named list of tidy DA tibbles, one per contrast. A single
#'   data.frame carrying a `contrast` column is also accepted and is split.
#' @param enrich_dfs Named list of tidy enrichment tibbles, one per contrast.
#'   Same split convenience as `volc_dfs`.
#' @param contrasts Character vector of contrast names to include and the
#'   order in which to draw them. Defaults to `names(volc_dfs)`.
#' @param subtitles Optional per-ring subtitles. Either a vector parallel to
#'   `contrasts`, or a vector named by contrast. `NULL` draws no subtitles.
#' @param nrow,ncol Outer layout dims; forwarded to `patchwork::wrap_plots()`.
#' @param tag_levels Panel-tag scheme; forwarded to `patchwork::plot_annotation()`.
#' @param guides Patchwork `guides` argument; default `"collect"` collects the
#'   shared NES legend.
#' @param panel_spacing Gutter between adjacent panels, in millimetres
#'   (default 1.5). Applied as half on each panel edge, so neighbours sit
#'   `panel_spacing` apart. Lower it to pack rings closer.
#' @param panel_margin Outer margin around the whole grid, in millimetres
#'   (default 2). Trims the dead frame around the assembled figure.
#' @param label_headroom Radial room reserved for pathway labels inside each
#'   panel (default 1.1, looser than the single-ring default of 0.5 so wide
#'   boxes stay enclosed when panels are packed tight). Forwarded to
#'   `volcano_ring()`; raise it if labels clip, lower it to enlarge the rings.
#' @param legend_position Placement of the collected NES legend: `"bottom"`
#'   (default, recovers the right-hand gap), `"right"`, or `"none"`.
#' @param legend_width Length of the NES colourbar long axis, in millimetres
#'   (default 26, tuned for the bottom bar). Sets the key width when the legend
#'   is horizontal, the key height when vertical; a side legend usually wants a
#'   larger value (~40).
#' @param ... Forwarded to each `volcano_ring()` call (e.g. `gene_col`,
#'   `padj_col`, `magnitude`, `theme`).
#' @return An S3 object `c("volcano_ring_grid", "list")` with elements
#'   `$plot` (patchwork) and `$data` (list of `list(volc, enrich)` pairs).
#' @export
volcano_ring_grid <- function(volc_dfs, enrich_dfs,
                              contrasts = NULL,
                              subtitles = NULL,
                              nrow = NULL,
                              ncol = NULL,
                              tag_levels = "A",
                              guides = "collect",
                              panel_spacing = 1.5,
                              panel_margin = 2,
                              label_headroom = 1.1,
                              legend_position = c("bottom", "right", "none"),
                              legend_width = 26,
                              ...) {
  legend_position <- match.arg(legend_position)
  validate_grid_spacing(panel_spacing, panel_margin, legend_width)
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

  panels <- lapply(seq_along(contrasts), function(i) {
    cn <- contrasts[[i]]
    sub <- if (is.null(subtitles)) {
      NULL
    } else if (!is.null(names(subtitles))) {
      subtitles[[cn]]
    } else {
      subtitles[[i]]
    }
    volcano_ring(volc_list[[cn]], enrich_list[[cn]],
      title = cn, subtitle = sub, label_headroom = label_headroom, ...
    )
  })
  names(panels) <- contrasts

  horizontal <- legend_position %in% c("bottom", "top")
  legend_theme <- ggplot2::theme(
    legend.position = legend_position,
    legend.direction = if (horizontal) "horizontal" else "vertical",
    legend.box.spacing = ggplot2::unit(panel_spacing, "mm"),
    legend.key.width = ggplot2::unit(if (horizontal) legend_width else 2, "mm"),
    legend.key.height = ggplot2::unit(if (horizontal) 2 else legend_width, "mm"),
    plot.margin = ggplot2::margin(
      panel_spacing / 2, panel_spacing / 2,
      panel_spacing / 2, panel_spacing / 2, "mm"
    )
  )
  legend_dir <- if (horizontal) "horizontal" else "vertical"

  plot <- patchwork::wrap_plots(panels, nrow = nrow, ncol = ncol, guides = guides) +
    patchwork::plot_annotation(
      tag_levels = tag_levels,
      theme = ggplot2::theme(
        plot.margin = ggplot2::margin(
          panel_margin, panel_margin, panel_margin, panel_margin, "mm"
        )
      )
    ) &
    legend_theme &
    ggplot2::guides(fill = ggplot2::guide_colorbar(direction = legend_dir))

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
