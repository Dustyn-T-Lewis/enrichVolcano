#' Build an enrichment ring ggplot for one contrast
#'
#' Uses `ggforce::geom_arc_bar` to display the top-N pathways as arc
#' segments. Magnitude maps to arc thickness, color to NES (or other),
#' terms split by up/down NES direction.
#'
#' @param enrich_result Tibble from `ev_collapse()` or `ev_enrich()`.
#' @param contrast Single contrast name.
#' @param max_terms Cap on terms shown (default 10; warn at 12, abort at 25).
#' @param order_by `"padj"` (default), `"size"`, `"NES"`, or `"alphabetic"`.
#' @param order_dir `"asc"` or `"desc"`.
#' @param magnitude What arc-thickness encodes: `"neg_log_padj"` (default),
#'   `"nes"`, or `"gene_count"`.
#' @param color Arc color encoding: `"nes"` (default), `"neg_log_padj"`,
#'   or `"direction"`.
#' @param split_by_direction Split up- and down-NES into two arc groups.
#' @param label_clean_fn Optional user fn applied to pathway names.
#' @param theme Output of `ev_theme()`.
#' @return ggplot object.
#' @export
ring_plot <- function(enrich_result, contrast,
                    max_terms = 10,
                    order_by = c("padj", "size", "NES", "alphabetic"),
                    order_dir = c("asc", "desc"),
                    magnitude = c("neg_log_padj", "nes", "gene_count"),
                    color = c("nes", "neg_log_padj", "direction"),
                    split_by_direction = TRUE,
                    label_clean_fn = NULL,
                    theme = ev_theme()) {
  order_by <- match.arg(order_by)
  order_dir <- match.arg(order_dir)
  magnitude <- match.arg(magnitude)
  color <- match.arg(color)

  if (max_terms > 25) {
    ev_abort("max_terms = {max_terms} exceeds legibility cap (25).",
             class = "ev_max_terms_exceeded")
  }
  if (max_terms > 12) {
    ev_warn("max_terms = {max_terms} > 12; ring labels may collide.",
            class = "ev_max_terms_warn")
  }

  sub <- enrich_result[enrich_result$contrast == contrast, , drop = FALSE]
  if ("dedup_kept" %in% colnames(sub)) {
    sub <- sub[sub$dedup_kept, , drop = FALSE]
  }

  if (nrow(sub) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0,
                          label = "No significant pathways", size = 5) +
        ggplot2::theme_void()
    )
  }

  sort_col <- switch(order_by,
    padj = "padj", size = "size",
    NES = "NES", alphabetic = "pathway"
  )
  sub <- sub[order(sub[[sort_col]],
                   decreasing = order_dir == "desc"), ]
  sub <- utils::head(sub, max_terms)
  sub$label <- if (!is.null(label_clean_fn)) {
    label_clean_fn(sub$pathway)
  } else {
    ev_clean_label(sub$pathway)
  }
  sub$mag <- switch(magnitude,
    neg_log_padj = -log10(sub$padj),
    nes          = sub$NES,
    gene_count   = sub$size
  )
  sub$col_val <- switch(color,
    nes          = sub$NES,
    neg_log_padj = -log10(sub$padj),
    direction    = ifelse(sub$direction == "up", 1, -1)
  )
  n <- nrow(sub)
  if (split_by_direction) {
    up <- sub[sub$direction == "up" & !is.na(sub$direction), ]
    dn <- sub[sub$direction == "down" & !is.na(sub$direction), ]
    if (nrow(up) > 0) {
      up$ang_start <- seq(0, pi * (nrow(up) - 1) / nrow(up),
                          length.out = nrow(up))
      up$ang_end <- up$ang_start + pi / nrow(up)
    }
    if (nrow(dn) > 0) {
      dn$ang_start <- seq(pi, pi + pi * (nrow(dn) - 1) / nrow(dn),
                          length.out = nrow(dn))
      dn$ang_end <- dn$ang_start + pi / nrow(dn)
    }
    sub <- dplyr::bind_rows(up, dn)
  } else {
    sub$ang_start <- seq(0, 2 * pi * (n - 1) / n, length.out = n)
    sub$ang_end <- sub$ang_start + 2 * pi / n
  }
  pal <- theme$palette
  ggplot2::ggplot(sub) +
    ggforce::geom_arc_bar(
      ggplot2::aes(x0 = 0, y0 = 0,
                   r0 = 1,
                   r = 1 + scales::rescale(.data$mag, to = c(0.1, 1)),
                   start = .data$ang_start,
                   end = .data$ang_end,
                   fill = .data$col_val)
    ) +
    ggplot2::scale_fill_gradient2(
      low = pal$nes_scale[1], mid = pal$nes_scale[3],
      high = pal$nes_scale[5], midpoint = 0
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::labs(title = contrast, fill = color)
}
