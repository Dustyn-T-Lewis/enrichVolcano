#' Build a volcano ggplot for one contrast
#'
#' @param data Tidy long tibble (post `pi_score`/`adjust_p`).
#' @param contrast Single contrast name.
#' @param p_method Which p-column to use on the y-axis: `"pi_eq2"`,
#'   `"pi_eq1"`, `"raw_p"`, or `"adj_p"`.
#' @param p_threshold Significance cutoff (default 0.05).
#' @param logfc_threshold Effect-size cutoff (default log2(1.5)).
#' @param label_mode Labeling rule: `"none"` (default), `"top_per_direction"`,
#'   `"top_total"`, `"all_significant"`, or `"explicit"`.
#' @param label_n Cap for the top-* modes (per direction for
#'   `"top_per_direction"`).
#' @param label_rank_by `"significance"` (default) or `"logfc"`.
#' @param label_genes Symbols/accessions to label when `label_mode =
#'   "explicit"`.
#' @param theme Output of `ev_theme()`.
#' @return ggplot object.
#' @export
ev_volcano <- function(data, contrast,
                       p_method = c("pi_eq2", "pi_eq1", "raw_p", "adj_p"),
                       p_threshold = 0.05,
                       logfc_threshold = log2(1.5),
                       label_mode = c("none", "top_per_direction", "top_total",
                                      "all_significant", "explicit"),
                       label_n = 10,
                       label_rank_by = c("significance", "logfc"),
                       label_genes = NULL,
                       theme = ev_theme()) {
  p_method <- match.arg(p_method)
  label_mode <- match.arg(label_mode)
  label_rank_by <- match.arg(label_rank_by)
  sub <- data[data$contrast == contrast, , drop = FALSE]
  if (nrow(sub) == 0) {
    ev_abort("No rows for contrast {.val {contrast}}.",
             class = "ev_unknown_contrast")
  }
  p_col <- switch(p_method,
    pi_eq2 = "pi_eq2", pi_eq1 = "pi_eq1",
    raw_p = "P.Value", adj_p = "adj.P.Val"
  )
  sub$ev_y <- -log10(sub[[p_col]])
  sub$ev_class <- dplyr::case_when(
    sub[[p_col]] < p_threshold & sub$logFC >=  logfc_threshold ~ "up",
    sub[[p_col]] < p_threshold & sub$logFC <= -logfc_threshold ~ "down",
    TRUE ~ "ns"
  )
  pal <- theme$palette
  top <- ev_select_labels(
    sub, mode = label_mode, n = label_n, rank_by = label_rank_by,
    genes = label_genes, p_col = p_col,
    p_threshold = p_threshold, logfc_threshold = logfc_threshold
  )
  ggplot2::ggplot(sub, ggplot2::aes(.data$logFC, .data$ev_y,
                                     colour = .data$ev_class)) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::scale_colour_manual(
      values = c(up = pal$up, down = pal$down, ns = pal$ns),
      breaks = c("up", "down", "ns")
    ) +
    ggrepel::geom_text_repel(
      data = top,
      ggplot2::aes(label = .data$label_text),
      max.overlaps = 50, force = 6, force_pull = 0.3,
      box.padding = 0.3, point.padding = 0.2, seed = 42
    ) +
    ggplot2::labs(
      x = "log2 fold change",
      y = paste0("-log10(", p_col, ")"),
      title = contrast,
      colour = NULL
    ) +
    theme$theme
}
