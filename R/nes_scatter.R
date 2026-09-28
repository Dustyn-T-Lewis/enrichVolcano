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
#'   `colour_by`. By significance, the fill is the theme's `up` colour for
#'   terms significant in `x` only, `down` for `y` only and `both` for both.
#' * Shading: same-direction quadrants take a tint of the `up` colour and
#'   opposite-direction quadrants a tint of the `down` colour.
#' * Corner labels: each quadrant's name and its count of significant terms,
#'   in the colour of its tint.
#' * Dashed line: \eqn{y = x}, or \eqn{y = -x} for a reversal.
#' * Labels: the most significant terms with at least `label_min_size` genes.
#'   The ten most significant are bold, and leader lines take the colour of
#'   their point.
#' * Subtitle: Spearman's \eqn{\rho} across all plotted terms, with a 95%
#'   Fieller interval when the correlation package is installed, its p-value, the
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
#' @param labels Term labels: a style, your own names or a function, as in
#'   [plot_volcano_ring()].
#' @param theme Output of [plot_theme()]; supplies the font and colours.
#'
#' @return A ggplot.
#' @references
#' Fieller EC, Hartley HO, Pearson ES (1957). Tests for rank correlation
#' coefficients. I. Biometrika 44(3-4):470-481. \doi{10.1093/biomet/44.3-4.470}
#'
#' Makowski D, Ben-Shachar MS, Patil I, Luedecke D (2020). Methods and
#' algorithms for correlation analysis in R. Journal of Open Source Software
#' 5(51):2306. \doi{10.21105/joss.02306}
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
                         labels = "short",
                         theme = plot_theme()) {
  check_enrichment(enrichment)
  check_labels(labels)
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
  share <- if (sum(counts) > 0) unname((counts[["top_right"]] + counts[["bottom_left"]]) / sum(counts)) else NA_real_
  if (comparison == "reversal") share <- 1 - share
  subtitle <- scatter_subtitle(
    score_correlation(wide$score_x, wide$score_y),
    share, c(concordance = "concordant", reversal = "reversed")[[comparison]],
    nrow(wide), sum(sig)
  )
  score_type <- enrichment@metadata$score_type

  draw_scatter(wide, sig, counts, comparison, colour_by, shape_by, label_min_size, label_n, labels, theme) +
    ggplot2::labs(
      x = sprintf("%s (%s)", score_type, x),
      y = sprintf("%s (%s)", score_type, y),
      subtitle = if (startsWith(subtitle, "rho")) bquote(rho ~ .(trimws(substring(subtitle, 4)))) else subtitle
    )
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
    # Adding 0 turns a rounded -0 into 0, so it prints without a minus sign.
    two <- function(x) sprintf("%.2f", round(x, 2) + 0)
    rho <- if (is.null(cor$ci)) {
      sprintf("rho = %s, %s", two(cor$rho), p_text)
    } else {
      sprintf("rho = %s [%s, %s], %s", two(cor$rho), two(cor$ci[1]), two(cor$ci[2]), p_text)
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
                         label_min_size, label_n, labels, theme) {
  pal <- theme$palette
  lim <- max(abs(c(wide$score_x, wide$score_y)), 1, na.rm = TRUE) * 1.5
  # Quadrants run top-right, top-left, bottom-left, bottom-right; the first and
  # third move the same way in both contrasts.
  tint <- rep(c(pal$up, pal$down), 2)
  if (!is.null(shape_by) && all(is.na(wide[[shape_by]]))) shape_by <- NULL
  sized <- !anyNA(wide$size)

  points <- wide[sig, , drop = FALSE]
  mapping <- list(x = quote(.data$score_x), y = quote(.data$score_y))
  fixed <- list(colour = "black", alpha = 0.95, stroke = 0.4)
  if (sized) mapping$size <- quote(.data$size) else fixed$size <- 3
  if (!is.null(colour_by)) mapping$fill <- rlang::expr(.data[[!!colour_by]]) else fixed$fill <- pal$both
  if (!is.null(shape_by)) mapping$shape <- rlang::expr(.data[[!!shape_by]]) else fixed$shape <- 21
  faint <- list(fill = pal$ns, colour = "grey55", alpha = 0.4, size = 1, stroke = 0.2)
  if (is.null(shape_by)) faint$shape <- 21

  # Every significant point goes to ggrepel, unlabelled ones with "", so labels
  # avoid all of them; point.size gives each its drawn size.
  size_range <- c(2, 8)
  rank <- order(pmin(points$sig_x, points$sig_y, na.rm = TRUE))
  eligible <- rank[is.na(points$size[rank]) | points$size[rank] >= label_min_size]
  shown <- utils::head(eligible, label_n)
  points$label <- character(nrow(points))
  points$label[shown] <- display_labels(points$term[shown], labels, 18)
  points$face <- rep("plain", nrow(points))
  points$face[utils::head(shown, 10)] <- "bold"
  points$point_size <- if (sized) scales::area_pal(size_range)(scales::rescale(points$size)) else rep(3, nrow(points))
  class_colours <- c(pal$up, pal$down, pal$both, pal$ns)
  segment <- if (identical(colour_by, "significance")) class_colours[points$significance] else "grey40"

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
      xmin = c(0, -lim, -lim, 0), xmax = c(lim, 0, 0, lim),
      ymin = c(0, 0, -lim, -lim), ymax = c(lim, lim, 0, 0),
      fill = tint, alpha = 0.12, colour = "grey70", linewidth = 0.2
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.2) +
    ggplot2::geom_abline(
      slope = if (comparison == "reversal") -1 else 1, intercept = 0,
      linetype = "dashed", linewidth = 0.3
    ) +
    do.call(ggplot2::geom_point, c(list(
      data = wide[!sig, , drop = FALSE],
      mapping = ggplot2::aes(!!!mapping[intersect(names(mapping), c("x", "y", "shape"))])
    ), faint)) +
    do.call(ggplot2::geom_point, c(list(data = points, mapping = ggplot2::aes(!!!mapping)), fixed)) +
    ggrepel::geom_text_repel(
      data = points,
      ggplot2::aes(x = .data$score_x, y = .data$score_y, label = .data$label, point.size = .data$point_size),
      fontface = points$face, segment.colour = segment, size = 2.6, lineheight = 0.85,
      bg.colour = "white", bg.r = 0.12, xlim = c(-lim, lim) * 0.93, ylim = c(-lim, lim) * 0.93,
      min.segment.length = 0, segment.size = 0.35, max.overlaps = Inf, box.padding = 0.4,
      point.padding = 0.2, force = 3, force_pull = 0.5, max.iter = 1e5, max.time = 3, seed = 42
    ) +
    ggplot2::annotate("label",
      x = c(lim, -lim, -lim, lim) * 0.98, y = c(lim, lim, -lim, -lim) * 0.98,
      label = paste0(quadrant_names[[comparison]], "  n = ", counts),
      hjust = c(1, 0, 0, 1), vjust = c(1, 1, 0, 0), colour = tint, fontface = "bold", size = 2.8,
      fill = scales::alpha("white", 0.92), linewidth = 0.3, label.padding = ggplot2::unit(1.5, "pt")
    ) +
    ggplot2::coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    ggplot2::theme_bw(base_size = theme$base_size, base_family = theme$base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold", colour = "grey30"),
      plot.subtitle = ggplot2::element_text(colour = "grey30"),
      legend.position = "bottom", legend.box = "vertical", legend.spacing.y = ggplot2::unit(0, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    )

  if (sized) {
    p <- p + ggplot2::scale_size_continuous(
      range = size_range, name = "Set size", guide = ggplot2::guide_legend(order = 3)
    )
  }
  if (!is.null(shape_by)) {
    n_shapes <- length(unique(stats::na.omit(wide[[shape_by]])))
    if (n_shapes > 5) {
      ev_abort("{.arg shape_by} has {n_shapes} values; point shapes support up to 5.",
        class = "enrichVolcano_param_error"
      )
    }
    p <- p + ggplot2::scale_shape_manual(
      values = c(21, 24, 22, 23, 25), na.value = 21, name = stringr::str_to_sentence(shape_by),
      guide = ggplot2::guide_legend(order = 2, override.aes = list(fill = "grey50", size = 3, alpha = 1))
    )
  }
  if (!is.null(colour_by)) p <- p + scatter_fill_scale(points[[colour_by]], colour_by, class_colours)
  if (!is.null(shape_by) && identical(colour_by, shape_by)) {
    p <- p + ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, override.aes = list(size = 3)),
      shape = ggplot2::guide_legend(order = 1)
    )
  }
  p
}

scatter_fill_scale <- function(values, name, class_colours) {
  if (is.numeric(values)) {
    return(ggplot2::scale_fill_viridis_c(name = stringr::str_to_sentence(name)))
  }
  guide <- ggplot2::guide_legend(order = 1, override.aes = list(shape = 21, size = 3))
  if (is.factor(values) && name == "significance") {
    colours <- stats::setNames(class_colours, levels(values))
    return(ggplot2::scale_fill_manual(values = colours, name = "Significant in", drop = TRUE, guide = guide))
  }
  name <- stringr::str_to_sentence(name)
  n <- length(unique(values))
  if (n <= 8) {
    okabe <- unname(grDevices::palette.colors(n + 1, palette = "Okabe-Ito")[-1])
    return(ggplot2::scale_fill_manual(values = okabe, name = name, guide = guide))
  }
  ggplot2::scale_fill_discrete(name = name, guide = guide)
}
