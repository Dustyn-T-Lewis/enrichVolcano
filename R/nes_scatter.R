#' Compare two contrasts term by term
#'
#' Plots each term's score in one contrast against its score in another.
#' Terms in the top-right and bottom-left quadrants move the same way in both
#' contrasts; terms in the other two move in opposite directions.
#'
#' `comparison` sets the question the figure asks:
#' * `"concordance"`: do two conditions regulate the same pathways, e.g.
#'   training in young against training in old? Same-direction quadrants are
#'   "concordant", the reference line is \eqn{y = x}, and the subtitle reports
#'   the share of significant terms that are concordant.
#' * `"reversal"`: does one condition undo another, e.g. aging against
#'   training? Opposite-direction quadrants are "reversed", same-direction
#'   ones "exacerbated", the reference line is \eqn{y = -x}, and the subtitle
#'   reports the share reversed.
#'
#' @section What is drawn:
#' * Points: every term scored in both contrasts. Non-significant terms are
#'   small and grey; significant ones are sized by gene-set size and filled by
#'   `colour_by`.
#' * Shading: same-direction quadrants are tinted, the others left white.
#' * Corner labels: each quadrant's name and its count of significant terms.
#' * Dashed line: \eqn{y = x}, or \eqn{y = -x} for a reversal.
#' * Labels: the most significant terms with at least `label_min_size` genes.
#' * Subtitle: Spearman's \eqn{\rho} across all plotted terms, with a 95%
#'   interval when the correlation package is installed, its p-value, the
#'   share of significant terms that are concordant (or reversed), and the
#'   term counts.
#'
#' A term is significant in a contrast when its `padj` is below `term_threshold`.
#' If the object has no adjusted p-values at all, nominal p-values are used
#' and a note says so.
#'
#' @param enrichment An enrichment object from [as_enrichment()] or
#'   [run_enrichment()] that holds both contrasts.
#' @param x,y Contrast names for the horizontal and vertical axes.
#' @param comparison `"concordance"` or `"reversal"`; see Description.
#' @param databases Collections to draw from, e.g. `c("Hallmark", "GO Slim")`;
#'   `NULL` (default) draws from all. Skipped, with a note, when the object has
#'   no database labels.
#' @param collapse Hide terms that [dedup_terms()] flagged redundant in one
#'   contrast and kept in neither.
#' @param term_threshold A term is significant in a contrast when its `padj`
#'   is below this.
#' @param colour_by Column that fills significant points: `"significance"`
#'   (significant in `x` only, `y` only, or both), any per-term column such as
#'   `"database"`, or `NULL` for one colour. Per-term columns are read from
#'   the `x` contrast.
#' @param shape_by Column mapped to point shape (up to five values), or `NULL`.
#' @param label_min_size Smallest gene set that gets a label.
#' @param label_n Most term labels drawn.
#' @param theme Output of [plot_theme()]; supplies the base font.
#'
#' @return A ggplot.
#' @export
#' @examples
#' ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
#'   package = "enrichVolcano"
#' )))
#' plot_scatter(ex, "Training_Young", "Training_Old")
#' plot_scatter(ex, "Aging", "Training_Old", comparison = "reversal")
plot_scatter <- function(enrichment, x, y,
                         comparison = c("concordance", "reversal"),
                         databases = NULL,
                         collapse = TRUE,
                         term_threshold = 0.05,
                         colour_by = "significance",
                         shape_by = "database",
                         label_min_size = 15,
                         label_n = 20,
                         theme = plot_theme()) {
  check_enrichment(enrichment)
  comparison <- rlang::arg_match(comparison)
  res <- enrichment@results
  if (identical(x, y)) {
    ev_abort("{.arg x} and {.arg y} must be different contrasts.", class = "enrichVolcano_input_error")
  }
  missing <- setdiff(c(x, y), res$contrast)
  if (length(missing) > 0) {
    ev_abort(
      "No contrast {.val {missing}} in this object; it holds {.val {unique(res$contrast)}}.",
      class = "enrichVolcano_input_error"
    )
  }

  sig_col <- "padj"
  if (all(is.na(res$padj))) {
    ev_inform("No adjusted p-values, so significance uses nominal p.", class = "enrichVolcano_nominal_p")
    sig_col <- "p"
  }
  wide <- pair_contrasts(filter_view(res, databases, collapse = FALSE), x, y)
  if (collapse) wide <- wide[!hidden_by_collapse(wide$dedup_status_x, wide$dedup_status_y), , drop = FALSE]
  wide$sig_x <- wide[[paste0(sig_col, "_x")]]
  wide$sig_y <- wide[[paste0(sig_col, "_y")]]
  wide$significance <- scatter_significance(wide$sig_x, wide$sig_y, term_threshold, x, y)
  for (arg in c("colour_by", "shape_by")) {
    col <- get(arg)
    if (!is.null(col) && !col %in% names(wide)) ev_abort_missing_column(wide, col, arg, "enrichment")
  }

  sig <- wide$significance != "NS"
  counts <- quadrant_counts(wide$score_x, wide$score_y, sig)
  share <- concordance_fraction(counts)
  if (comparison == "reversal") share <- 1 - share
  subtitle <- scatter_subtitle(
    score_correlation(wide$score_x, wide$score_y),
    share, c(concordance = "concordant", reversal = "reversed")[[comparison]],
    nrow(wide), sum(sig)
  )
  score_type <- enrichment@metadata$score_type

  draw_scatter(wide, sig, counts, comparison, colour_by, shape_by, label_min_size, label_n, theme) +
    ggplot2::labs(
      x = sprintf("%s (%s)", score_type, x),
      y = sprintf("%s (%s)", score_type, y),
      subtitle = subtitle_math(subtitle)
    )
}

subtitle_math <- function(text) {
  if (!startsWith(text, "rho")) {
    return(text)
  }
  bquote(rho ~ .(trimws(substring(text, 4))))
}

pair_contrasts <- function(res, x, y) {
  a <- res[res$contrast == x, , drop = FALSE]
  b <- res[res$contrast == y, , drop = FALSE]
  key <- function(d) paste(d$database, d$term, sep = "\r")
  shared <- intersect(key(a), key(b))
  if (length(shared) == 0) {
    ev_abort("Contrasts {.val {x}} and {.val {y}} share no terms to compare.",
      class = "enrichVolcano_input_error"
    )
  }
  n_single <- length(union(key(a), key(b))) - length(shared)
  if (n_single > 0) {
    ev_inform("{n_single} term{?s} appear{?s/} in only one contrast and {?is/are} left out.",
      class = "enrichVolcano_unpaired_terms"
    )
  }
  a <- a[match(shared, key(a)), , drop = FALSE]
  b <- b[match(shared, key(b)), , drop = FALSE]
  per_contrast <- c("contrast", "score", "p", "padj", "direction", "leading_edge", "dedup_status")
  wide <- a[setdiff(names(a), per_contrast)]
  wide$size <- ifelse(is.na(a$size), b$size, a$size)
  for (col in c("score", "p", "padj", "dedup_status")) {
    wide[[paste0(col, "_x")]] <- if (col %in% names(a)) a[[col]] else NA
    wide[[paste0(col, "_y")]] <- if (col %in% names(b)) b[[col]] else NA
  }
  rownames(wide) <- NULL
  wide
}

hidden_by_collapse <- function(status_x, status_y) {
  redundant <- status_x %in% "redundant" | status_y %in% "redundant"
  kept <- status_x %in% "kept" | status_y %in% "kept"
  redundant & !kept
}

scatter_significance <- function(padj_x, padj_y, term_threshold, x_name, y_name) {
  in_x <- !is.na(padj_x) & padj_x < term_threshold
  in_y <- !is.na(padj_y) & padj_y < term_threshold
  levels <- c(paste(x_name, "only"), paste(y_name, "only"), "Both", "NS")
  class <- ifelse(in_x & in_y, "Both", ifelse(in_x, levels[1], ifelse(in_y, levels[2], "NS")))
  factor(class, levels = levels)
}

quadrant_counts <- function(score_x, score_y, sig) {
  c(
    top_right = sum(sig & score_x > 0 & score_y > 0),
    top_left = sum(sig & score_x < 0 & score_y > 0),
    bottom_left = sum(sig & score_x < 0 & score_y < 0),
    bottom_right = sum(sig & score_x > 0 & score_y < 0)
  )
}

concordance_fraction <- function(counts) {
  if (sum(counts) == 0) {
    return(NA_real_)
  }
  unname((counts[["top_right"]] + counts[["bottom_left"]]) / sum(counts))
}

score_correlation <- function(x, y, with_ci = requireNamespace("correlation", quietly = TRUE)) {
  if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NULL)
  }
  if (with_ci) {
    r <- correlation::cor_test(data.frame(x = x, y = y), "x", "y", method = "spearman")
    return(list(rho = r$rho, ci = c(r$CI_low, r$CI_high), p = r$p))
  }
  ev_inform("Install the correlation package to add a confidence interval for rho.",
    class = "enrichVolcano_no_ci"
  )
  r <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  list(rho = unname(r$estimate), ci = NULL, p = r$p.value)
}

scatter_subtitle <- function(cor, share, share_label, n, n_sig) {
  rho <- NULL
  if (!is.null(cor)) {
    p_text <- if (cor$p < 0.001) "p < 0.001" else sprintf("p = %.2f", cor$p)
    rho <- if (is.null(cor$ci)) {
      sprintf("rho = %.2f, %s", cor$rho, p_text)
    } else {
      sprintf("rho = %.2f [%.2f, %.2f], %s", cor$rho, cor$ci[1], cor$ci[2], p_text)
    }
  }
  paste(c(
    rho,
    if (!is.na(share)) sprintf("%.0f%% %s", 100 * share, share_label),
    sprintf("%d terms, %d significant", n, n_sig)
  ), collapse = " | ")
}

quadrant_names <- list(
  concordance = c("Concordant up", "Discordant", "Concordant down", "Discordant"),
  reversal = c("Exacerbated", "Reversed", "Exacerbated", "Reversed")
)

draw_scatter <- function(wide, sig, counts, comparison, colour_by, shape_by,
                         label_min_size, label_n, theme) {
  lim <- max(abs(c(wide$score_x, wide$score_y)), 1, na.rm = TRUE) * 1.35
  quad <- data.frame(
    xmin = c(0, -lim, -lim, 0), xmax = c(lim, 0, 0, lim),
    ymin = c(0, 0, -lim, -lim), ymax = c(lim, lim, 0, 0),
    fill = c("grey93", "white", "grey93", "white")
  )
  corner <- data.frame(
    x = c(lim, -lim, -lim, lim) * 0.97, y = c(lim, lim, -lim, -lim) * 0.97,
    hjust = c(1, 0, 0, 1), vjust = c(1, 1, 0, 0),
    label = paste0(quadrant_names[[comparison]], "\nn = ", counts)
  )
  if (!is.null(shape_by) && all(is.na(wide[[shape_by]]))) shape_by <- NULL
  sized <- !anyNA(wide$size)

  points <- wide[sig, , drop = FALSE]
  mapping <- list(x = quote(.data$score_x), y = quote(.data$score_y))
  fixed <- list(colour = "grey20", alpha = 0.85, stroke = 0.3)
  if (sized) mapping$size <- quote(.data$size) else fixed$size <- 2.5
  if (!is.null(colour_by)) mapping$fill <- rlang::expr(.data[[!!colour_by]]) else fixed$fill <- "#0072B2"
  if (!is.null(shape_by)) mapping$shape <- rlang::expr(.data[[!!shape_by]]) else fixed$shape <- 21

  labelled <- points[is.na(points$size) | points$size >= label_min_size, , drop = FALSE]
  labelled <- labelled[order(pmin(labelled$sig_x, labelled$sig_y, na.rm = TRUE)), , drop = FALSE]
  labelled <- utils::head(labelled, label_n)
  labelled$label <- clean_label(labelled$term, width = 20)

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
      xmin = quad$xmin, xmax = quad$xmax, ymin = quad$ymin, ymax = quad$ymax, fill = quad$fill
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.2) +
    ggplot2::geom_abline(
      slope = if (comparison == "reversal") -1 else 1, intercept = 0,
      linetype = "dashed", linewidth = 0.3
    ) +
    ggplot2::geom_point(
      data = wide[!sig, , drop = FALSE],
      ggplot2::aes(x = .data$score_x, y = .data$score_y),
      shape = 21, size = 1, fill = "grey80", colour = "grey60", alpha = 0.5, stroke = 0.2
    ) +
    do.call(ggplot2::geom_point, c(list(data = points, mapping = ggplot2::aes(!!!mapping)), fixed)) +
    ggrepel::geom_label_repel(
      data = labelled, ggplot2::aes(x = .data$score_x, y = .data$score_y, label = .data$label),
      size = 2.6, fill = "white", label.size = 0.2, lineheight = 0.85,
      min.segment.length = 0, max.overlaps = Inf, box.padding = 0.4, seed = 42
    ) +
    ggplot2::annotate("label",
      x = corner$x, y = corner$y, label = corner$label,
      hjust = corner$hjust, vjust = corner$vjust, fontface = "bold", size = 2.8, lineheight = 0.9,
      fill = "white", linewidth = 0
    ) +
    ggplot2::coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    ggplot2::theme_classic(base_size = theme$base_size, base_family = theme$base_family) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

  if (sized) {
    p <- p + ggplot2::scale_size_continuous(
      range = c(1.5, 5), name = "Set size", guide = ggplot2::guide_legend(order = 3)
    )
  }
  if (!is.null(shape_by)) {
    n_shapes <- length(unique(stats::na.omit(points[[shape_by]])))
    if (n_shapes > 5) {
      ev_abort("{.arg shape_by} has {n_shapes} values; point shapes support up to 5.",
        class = "enrichVolcano_param_error"
      )
    }
    p <- p + ggplot2::scale_shape_manual(
      values = c(21, 24, 22, 23, 25), na.value = 21, name = shape_by,
      guide = ggplot2::guide_legend(order = 2)
    )
  }
  if (!is.null(colour_by)) p <- p + scatter_fill_scale(points[[colour_by]], colour_by)
  p
}

scatter_fill_scale <- function(values, name) {
  if (is.numeric(values)) {
    return(ggplot2::scale_fill_viridis_c(name = name))
  }
  guide <- ggplot2::guide_legend(order = 1, override.aes = list(shape = 21, size = 3))
  if (is.factor(values) && name == "significance") {
    okabe <- unname(grDevices::palette.colors(palette = "Okabe-Ito"))
    colours <- stats::setNames(c(okabe[2:4], "grey80"), levels(values))
    return(ggplot2::scale_fill_manual(values = colours, name = "Significant in", drop = TRUE, guide = guide))
  }
  n <- length(unique(values))
  if (n <= 8) {
    okabe <- unname(grDevices::palette.colors(n + 1, palette = "Okabe-Ito")[-1])
    return(ggplot2::scale_fill_manual(values = okabe, name = name, guide = guide))
  }
  ggplot2::scale_fill_discrete(name = name, guide = guide)
}
