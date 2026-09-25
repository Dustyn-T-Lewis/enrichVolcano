# Volcano-in-ring composite: a single Cartesian ggplot with the volcano at the
# centre and enrichment pathways wrapped around it as an outer ring. Arcs are
# drawn with ggforce::geom_arc_bar() in a coord_fixed() panel.

ev_ring_r_inner <- 4.4
ev_ring_r_outer <- 4.8

#' @keywords internal
#' @noRd
ev_ring_geometry <- function(enrich_df, term_col, padj_col, nes_col,
                             magnitude_col, genes_list, order_by = "padj",
                             gap_intra = 3, gap_split = 8,
                             arc_r0 = ev_ring_r_outer,
                             min_height = 0.05, max_height = 1.6, proportional = FALSE,
                             labels = NULL) {
  if (nrow(enrich_df) == 0) {
    return(enrich_df)
  }

  ring <- enrich_df
  ring <- ring[!is.na(ring[[term_col]]), , drop = FALSE]
  ring[[padj_col]][is.na(ring[[padj_col]])] <- 1
  ring[[nes_col]][is.na(ring[[nes_col]])] <- 0
  if (nrow(ring) == 0) {
    return(ring)
  }
  ring$.ev_clean_label <- display_labels(ring[[term_col]], labels, 15)
  ring$.ev_genes <- genes_list
  ring$.ev_is_up <- ring[[nes_col]] > 0
  ring$.ev_magnitude <- magnitude_col

  up <- ring[ring$.ev_is_up, , drop = FALSE]
  dn <- ring[!ring$.ev_is_up, , drop = FALSE]
  ord_key <- function(df) {
    switch(order_by,
      padj = df[[padj_col]],
      score = -abs(df[[nes_col]])
    )
  }
  up <- up[order(ord_key(up)), , drop = FALSE]
  dn <- dn[order(ord_key(dn)), , drop = FALSE]
  n_up <- nrow(up)
  n_dn <- nrow(dn)

  scale_h <- function(df) {
    if (nrow(df) == 0) {
      return(numeric(0))
    }
    if (proportional) {
      return(arc_r0 + max_height * df$.ev_magnitude)
    }
    if (nrow(df) == 1) {
      return(arc_r0 + (min_height + max_height) / 2)
    }
    mag <- df$.ev_magnitude
    rng <- range(mag)
    if (diff(rng) <= 0) {
      return(rep(arc_r0 + (min_height + max_height) / 2, nrow(df)))
    }
    scaled <- (mag - rng[1]) / diff(rng)
    arc_r0 + min_height + (max_height - min_height) * sqrt(scaled)
  }
  if (n_up > 0) up$arc_r1_var <- scale_h(up)
  if (n_dn > 0) dn$arc_r1_var <- scale_h(dn)

  # Up terms fill the right half from the top, down terms the left half from the top.
  place <- function(df, centre, span, from_top) {
    n <- nrow(df)
    if (n == 0) {
      return(df)
    }
    arc_w <- (span - (n - 1) * gap_intra) / n
    slot <- if (from_top) seq_len(n) - 1 else n - seq_len(n)
    df$start_deg <- centre - span / 2 + slot * (arc_w + gap_intra)
    df$end_deg <- df$start_deg + arc_w
    df
  }
  if (n_up > 0 && n_dn > 0) {
    arc_w <- (360 - 2 * gap_split - (n_up - 1) * gap_intra - (n_dn - 1) * gap_intra) / (n_up + n_dn)
    up_span <- n_up * arc_w + (n_up - 1) * gap_intra
    dn_span <- n_dn * arc_w + (n_dn - 1) * gap_intra
  } else {
    half <- function(n) if (n == 1) 30 else 180 - gap_split
    up_span <- half(n_up)
    dn_span <- half(n_dn)
  }
  up <- place(up, 90, up_span, from_top = TRUE)
  dn <- place(dn, 270, dn_span, from_top = FALSE)

  ring <- rbind(up, dn)
  ring$.ev_mid_deg <- (ring$start_deg + ring$end_deg) / 2
  ring$.ev_start_rad <- ring$start_deg * pi / 180
  ring$.ev_end_rad <- ring$end_deg * pi / 180
  ring$.ev_mid_rad <- ring$.ev_mid_deg * pi / 180
  ring$.ev_term_idx <- seq_len(nrow(ring))
  ring
}

#' @keywords internal
#' @noRd
ev_tick_data <- function(ring_data, da, gene_col, logfc_col,
                         tick_r0 = ev_ring_r_inner, tick_r1 = ev_ring_r_outer) {
  if (nrow(ring_data) == 0) {
    return(data.frame())
  }
  gene_lfc <- da[!is.na(da[[logfc_col]]),
    c(gene_col, logfc_col),
    drop = FALSE
  ]
  names(gene_lfc) <- c("gene", "logFC")
  gene_lfc <- gene_lfc[!duplicated(gene_lfc$gene), , drop = FALSE]
  pad <- 0.5 * pi / 180

  out <- lapply(seq_len(nrow(ring_data)), function(i) {
    row <- ring_data[i, ]
    genes <- intersect(row$.ev_genes[[1]], gene_lfc$gene)
    if (length(genes) == 0) {
      return(NULL)
    }
    a0 <- row$.ev_start_rad + pad
    a1 <- row$.ev_end_rad - pad
    if (a1 <= a0) a1 <- a0 + pad
    ang <- seq(a0, a1, length.out = length(genes))
    lfc <- gene_lfc$logFC[match(genes, gene_lfc$gene)]
    data.frame(
      gene = genes,
      logFC = lfc,
      direction = ifelse(lfc > 0, "up", "down"),
      x0 = tick_r0 * sin(ang), y0 = tick_r0 * cos(ang),
      x1 = tick_r1 * sin(ang), y1 = tick_r1 * cos(ang),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), out)
  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}


#' Volcano-in-ring composite for one contrast
#'
#' Draws a differential-abundance volcano embedded in a ring of enrichment
#' terms. Score-coloured arcs sit around the volcano; tick lines drop from the
#' arcs to each term's leading-edge genes inside the volcano.
#'
#' @section Which terms ring the volcano:
#' From the chosen contrast, terms in `databases` are kept, redundant ones are
#' hidden when `collapse = TRUE` (see [dedup_terms()]), and the `n_terms` terms with
#' the smallest `padj` below `term_threshold` are drawn, whatever their
#' direction, so the up and down halves need not match. A term found in two
#' collections is drawn once. `terms` overrides all of this with a hand-picked
#' set.
#'
#' @param da DA results from [as_da()]; the rows for `contrast` are drawn.
#'   Points are called significant on `padj`, or on `p` when `padj` is empty.
#'   To colour by another statistic, such as a pi-value, pass it to [as_da()]
#'   as `padj`.
#' @param enrichment An enrichment object from [as_enrichment()] or [run_enrichment()].
#' @param contrast Which contrast of `enrichment` to draw. Needed only when it
#'   holds more than one.
#' @param databases Collections to draw from, matched against the `database`
#'   column, e.g. `c("Hallmark", "GO Slim")`. `NULL` (default) draws from all.
#'   Skipped, with a note, when the object has no database labels.
#' @param collapse Hide terms flagged `"redundant"` by [dedup_terms()].
#' @param term_threshold Terms need `padj` below this to be drawn.
#' @param n_terms Most unique terms drawn, counted across both directions.
#' @param terms Optional character vector of exact term names to draw instead.
#' @param labels Your own display names, as a character vector named by term,
#'   such as `c(HALLMARK_OXIDATIVE_PHOSPHORYLATION = "OXPHOS")`. Terms not
#'   named keep [clean_label()]. A name without `\n` is wrapped like the
#'   others.
#' @param p_threshold Significance cutoff for volcano points.
#' @param logfc_threshold Effect-size cutoff; a point is called up/down only
#'   when `abs(logFC) >= logfc_threshold` as well as significant.
#' @param title,subtitle,tag Title, subtitle and panel tag.
#' @param volcano_radius Radius of the volcano in plot units.
#' @param x_scale Horizontal compression of the point cloud (default 1). Values
#'   below 1 pull points toward the fold-change axis so the widest points clear
#'   the enrichment ring; the up/down axis annotations are unaffected.
#' @param y_scale Vertical compression of the point cloud (default 1), anchored
#'   at the fold-change axis. Values below 1 lower the tallest points so they
#'   clear the `-log10 p` label at the top of the volcano.
#' @param ring_radius Inner radius of the enrichment ring (default 4.8), where
#'   the leading-edge tick band begins. Raise it to widen the central breathing
#'   gap around the volcano, lower it to close it. Keep it above
#'   `volcano_radius * 0.92` so the point cloud clears the ring.
#' @param ring_thickness Radial width of the tick band between `ring_radius` and
#'   the foot of the coloured arcs (default 0.55). This is the length of the
#'   leading-edge ticks; widen it to make ticks easier to read.
#' @param tick_width Line width of the leading-edge ticks (default 0.3).
#' @param label_headroom Extra radial room (data units, default 0.5) reserved
#'   beyond the outermost pathway label so its box stays enclosed within the
#'   square panel rather than clipping or spilling into a neighbour. Raise it
#'   when wide label boxes are clipped; lower it to pack the ring tighter.
#' @param disc_colour Optional fill for a tinted central disc.
#' @param score_limits Length-2 numeric or `NULL`; limits of the arc fill scale.
#'   `NULL` uses `c(-3, 3)` for NES, and otherwise spans the largest absolute
#'   significant score in the whole object, so every contrast shares a scale.
#' @param magnitude What arc height encodes: `"neg_log_padj"` or `"size"`.
#'   `NULL` picks `"neg_log_padj"` for NES and `"size"` otherwise, so fill and
#'   height never repeat the same number for fry or camera results; without
#'   set sizes it falls back to `"neg_log_padj"`.
#' @param arc_order Angular order of arcs within each up/down half:
#'   `"padj"` (default, lowest FDR first) or `"score"` (strongest absolute score
#'   first). The up/down split itself is always by direction.
#' @param arc_height_range Length-2 numeric `c(min, max)` for the shortest
#'   and tallest arc; widen it to exaggerate the magnitude encoding.
#' @param show_counts Draw the up/down significant-point count badges.
#' @param point_size,point_alpha Volcano point size (default 1.1) and opacity.
#' @param label_size Text size of the term labels and the point labels.
#' @param label_gap Radial gap between each arc's outer top and its own label
#'   (default 0.6). Anchoring per-arc keeps every leader line the same short
#'   length regardless of arc height; the label box grows outward from this
#'   point so it never overlaps the arc. Widen it to lengthen all leaders.
#' @param count_size Text size of the up/down count badges (default 2.4).
#' @param count_x_mult,count_y_mult Badge position as a fraction of the volcano
#'   radius (default 0.7).
#' @param axis_size Text size of the `up`/`down`/`log2 FC`/`-log10 p` axis
#'   annotations (default 2.2).
#' @param label_mode Which points get a gene label: `"none"`,
#'   `"top_per_direction"` (`label_n` up and `label_n` down),
#'   `"by_significance"` (`label_n` in total) or `"by_genes"` (the proteins in
#'   `label_genes`).
#' @param label_n How many points the two top modes label.
#' @param label_rank_by Rank points for the top modes by `"significance"`
#'   (smallest `p`) or `"logfc"` (largest absolute fold change).
#' @param label_genes Gene symbols or accessions to label with `"by_genes"`.
#' @param theme Output of [plot_theme()].
#' @return A ggplot.
#' @export
#' @examples
#' da <- read_example("yvo")$da
#' ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
#'   package = "enrichVolcano"
#' )))
#'
#' plot_volcano_ring(da, ex, contrast = "Training_Young", title = "Training_Young")
plot_volcano_ring <- function(da, enrichment,
                              contrast = NULL,
                              databases = NULL,
                              collapse = TRUE,
                              term_threshold = 0.05,
                              n_terms = 12,
                              terms = NULL,
                              labels = NULL,
                              p_threshold = 0.05,
                              logfc_threshold = 0,
                              title = NULL,
                              subtitle = NULL,
                              tag = NULL,
                              volcano_radius = 4.0,
                              x_scale = 1,
                              y_scale = 1,
                              ring_radius = 4.8,
                              ring_thickness = 0.55,
                              tick_width = 0.3,
                              label_headroom = 0.5,
                              disc_colour = NULL,
                              score_limits = NULL,
                              magnitude = NULL,
                              arc_order = c("padj", "score"),
                              arc_height_range = c(0.4, 1.6),
                              show_counts = TRUE,
                              point_size = 1.1,
                              point_alpha = 0.85,
                              label_size = 2.8,
                              label_gap = 0.6,
                              count_size = 2.4,
                              count_x_mult = 0.7,
                              count_y_mult = 0.7,
                              axis_size = 2.2,
                              label_mode = c(
                                "none", "top_per_direction",
                                "by_significance", "by_genes"
                              ),
                              label_n = 5,
                              label_rank_by = c("significance", "logfc"),
                              label_genes = NULL,
                              theme = plot_theme()) {
  do.call(draw_ring, as.list(environment()))
}

# The body of plot_volcano_ring(), shared with plot_bias_ring(). It takes every
# plot_volcano_ring() argument; `normalise` scales fill and arc height to the
# strongest drawn term.
draw_ring <- function(da, enrichment, contrast, databases, collapse, term_threshold, n_terms, terms, labels,
                      p_threshold, logfc_threshold, title, subtitle, tag, volcano_radius, x_scale,
                      y_scale, ring_radius, ring_thickness, tick_width, label_headroom, disc_colour,
                      score_limits, magnitude, arc_order, arc_height_range, show_counts, point_size,
                      point_alpha, label_size, label_gap, count_size, count_x_mult, count_y_mult,
                      axis_size, label_mode, label_n, label_rank_by, label_genes, theme,
                      term_labels = TRUE, normalise = FALSE) {
  arc_order <- match.arg(arc_order, c("padj", "score"))
  label_mode <- match.arg(label_mode, c("none", "top_per_direction", "by_significance", "by_genes"))
  label_rank_by <- match.arg(label_rank_by, c("significance", "logfc"))
  ev_assert_colour(disc_colour)
  validate_ring_geometry(ring_radius, volcano_radius, arc_height_range, label_headroom)

  if (!inherits(da, "enrichVolcano_da")) {
    ev_abort(
      "{.arg da} must be {.fn as_da} output, not {.cls {class(da)[1]}}.",
      class = "enrichVolcano_input_error"
    )
  }
  check_enrichment(enrichment)
  check_labels(labels)
  require_columns(da, c("gene", "logFC", "p", "padj"))
  validate_da(da)
  da <- contrast_rows(da, contrast)
  if (all(is.na(da$p)) && all(is.na(da$padj))) {
    ev_abort("The DA results have no p-values to draw a volcano from.", class = "enrichVolcano_data_error")
  }
  pval_col <- if (all(is.na(da$p))) "padj" else "p"
  score_type <- enrichment@metadata$score_type
  enrich_df <- contrast_rows(enrichment@results, contrast)
  if (is.null(terms)) enrich_df <- filter_view(enrich_df, databases, collapse)
  enrich_df <- ring_terms(enrich_df, term_threshold, n_terms, terms)
  magnitude <- default_magnitude(magnitude, score_type, enrich_df$size)
  if (nrow(enrich_df) == 0) {
    ev_inform("No terms pass the selection, so the ring is empty.", class = "enrichVolcano_empty_ring")
  }
  has_padj <- !all(is.na(da$padj))

  pal <- theme$palette
  score_limits <- score_limits %||% theme$score_limits %||%
    default_score_limits(score_type, enrichment@results, term_threshold)
  if (normalise) {
    subtitle <- subtitle %||% sprintf(
      "%d up, %d down of %d significant terms",
      sum(enrich_df$score > 0), sum(enrich_df$score < 0), nrow(enrich_df)
    )
    if (nrow(enrich_df) > 0) enrich_df$score <- enrich_df$score / max(abs(enrich_df$score))
    score_type <- paste(score_type, "/ max")
    score_limits <- c(-1, 1)
  }
  vr <- volcano_radius * 0.92
  ring_r0 <- ring_radius
  ring_r1 <- ring_radius + ring_thickness

  v <- da[!is.na(da$logFC) & !is.na(da[[pval_col]]), ,
    drop = FALSE
  ]
  v$.ev_neg_log10p <- -log10(v[[pval_col]])
  v <- v[is.finite(v$.ev_neg_log10p), , drop = FALSE]

  sig <- (if (has_padj) v$padj else v[[pval_col]]) < p_threshold
  lfc <- v$logFC
  v$.ev_direction <- ifelse(sig & lfc >= logfc_threshold, "up",
    ifelse(sig & lfc <= -logfc_threshold, "down", "ns")
  )
  n_up <- sum(v$.ev_direction == "up")
  n_down <- sum(v$.ev_direction == "down")

  x_max <- max(abs(lfc), na.rm = TRUE)
  y_max <- max(v$.ev_neg_log10p, na.rm = TRUE)
  if (!is.finite(x_max) || x_max == 0) x_max <- 1
  if (!is.finite(y_max) || y_max == 0) y_max <- 1
  v$.ev_x_plot <- lfc / x_max * vr * x_scale
  v$.ev_y_plot <- v$.ev_neg_log10p / y_max * 2 * vr * y_scale - vr
  v_ns <- v[v$.ev_direction == "ns", , drop = FALSE]
  v_sig <- v[v$.ev_direction != "ns", , drop = FALSE]

  v_labels <- v
  v_labels$x_plot <- v$.ev_x_plot
  v_labels$y_plot <- v$.ev_y_plot
  lab_pts <- ev_select_labels(
    v_labels,
    mode = label_mode, n = label_n, rank_by = label_rank_by,
    genes = label_genes, p_col = pval_col,
    p_threshold = p_threshold, logfc_threshold = logfc_threshold
  )

  mag_vec <- if (normalise) {
    abs(enrich_df$score)
  } else {
    switch(magnitude,
      neg_log_padj = -log10(pmax(enrich_df$padj, .Machine$double.xmin)),
      size         = enrich_df$size
    )
  }
  ring <- ev_ring_geometry(enrich_df,
    term_col = "term", padj_col = "padj",
    nes_col = "score", magnitude_col = mag_vec,
    genes_list = enrich_df$leading_edge, order_by = arc_order, arc_r0 = ring_r1,
    min_height = arc_height_range[1], max_height = arc_height_range[2], proportional = normalise,
    labels = labels
  )
  ticks <- ev_tick_data(ring, v,
    gene_col = "gene", logfc_col = "logFC",
    tick_r0 = ring_r0, tick_r1 = ring_r1
  )

  p <- ggplot2::ggplot()

  if (!is.null(disc_colour)) {
    disc <- data.frame(
      x = ring_r0 * cos(seq(0, 2 * pi, length.out = 200)),
      y = ring_r0 * sin(seq(0, 2 * pi, length.out = 200))
    )
    p <- p + ggplot2::geom_polygon(
      data = disc, ggplot2::aes(.data$x, .data$y),
      fill = disc_colour, alpha = 0.12, colour = NA, inherit.aes = FALSE
    )
  }

  if (nrow(ring) > 0) {
    p <- p + ggforce::geom_arc_bar(
      data = ring,
      ggplot2::aes(
        x0 = 0, y0 = 0, r0 = ring_r0, r = ring_r1,
        start = .data$.ev_start_rad, end = .data$.ev_end_rad
      ),
      fill = "grey93", colour = "grey78", linewidth = 0.15,
      inherit.aes = FALSE
    )
  }

  p <- p +
    ggplot2::geom_point(
      data = v_ns, ggplot2::aes(.data$.ev_x_plot, .data$.ev_y_plot),
      colour = pal$ns, size = point_size * 0.7, alpha = point_alpha * 0.4
    ) +
    ggplot2::geom_point(
      data = v_sig,
      ggplot2::aes(.data$.ev_x_plot, .data$.ev_y_plot, colour = .data$.ev_direction),
      size = point_size, alpha = point_alpha
    ) +
    ggplot2::scale_colour_manual(
      values = c(up = pal$up, down = pal$down), guide = "none"
    ) +
    ggplot2::annotate("segment",
      x = 0, xend = 0, y = -vr, yend = vr * 0.96,
      linewidth = 0.3, linetype = "dashed", colour = "grey50",
      arrow = ggplot2::arrow(
        length = ggplot2::unit(1.2, "mm"),
        type = "closed"
      )
    ) +
    ggplot2::annotate("segment",
      x = -vr * 0.42, xend = vr * 0.42,
      y = -vr, yend = -vr, linewidth = 0.3,
      linetype = "dashed", colour = "grey50",
      arrow = ggplot2::arrow(
        ends = "both",
        length = ggplot2::unit(1.2, "mm"),
        type = "closed"
      )
    ) +
    ggplot2::annotate("text",
      x = 0, y = vr * 1.05,
      label = "-log10 p", size = axis_size, colour = "grey40",
      parse = FALSE
    ) +
    ggplot2::annotate("text",
      x = vr * 0.45, y = -vr, label = "up",
      size = axis_size, colour = pal$up,
      fontface = "bold.italic", hjust = 0
    ) +
    ggplot2::annotate("text",
      x = -vr * 0.45, y = -vr, label = "down",
      size = axis_size, colour = pal$down,
      fontface = "bold.italic", hjust = 1
    ) +
    ggplot2::annotate("text",
      x = 0, y = -vr - 0.35, label = "log2 FC",
      size = axis_size, colour = "grey40", fontface = "bold.italic"
    )

  if (show_counts) {
    count_badge <- function(p, x, fill) {
      p +
        ggplot2::annotate("label",
          x = x, y = vr * count_y_mult, label = if (x > 0) n_up else n_down,
          size = count_size, colour = "grey20", fill = fill, fontface = "bold",
          linewidth = 0.3, label.r = ggplot2::unit(2, "pt"),
          label.padding = ggplot2::unit(2, "pt")
        ) +
        ggplot2::annotate("text",
          x = x, y = vr * count_y_mult, label = if (x > 0) n_up else n_down,
          size = count_size, colour = "white", fontface = "bold"
        )
    }
    p <- p |>
      count_badge(vr * count_x_mult, pal$up) |>
      count_badge(-vr * count_x_mult, pal$down)
  }

  max_r <- 5.6
  x_half <- NULL
  y_half <- NULL
  if (nrow(ring) > 0) {
    if (nrow(ticks) > 0) {
      p <- p + ggplot2::geom_segment(
        data = ticks,
        ggplot2::aes(
          x = .data$x0, y = .data$y0, xend = .data$x1,
          yend = .data$y1, colour = .data$direction
        ),
        linewidth = tick_width, alpha = 0.8, inherit.aes = FALSE
      )
    }
    p <- p +
      ggforce::geom_arc_bar(
        data = ring,
        ggplot2::aes(
          x0 = 0, y0 = 0, r0 = ring_r1, r = .data$arc_r1_var,
          start = .data$.ev_start_rad, end = .data$.ev_end_rad,
          fill = .data$score
        ),
        colour = "grey40", linewidth = 0.2, inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_gradientn(
        colours = pal$nes_scale,
        values = scales::rescale(pal$nes_values),
        limits = score_limits, oob = scales::squish, name = score_type
      )
  }
  if (nrow(ring) > 0 && !term_labels) max_r <- max(ring$arc_r1_var) + label_headroom
  if (nrow(ring) > 0 && term_labels) {
    lbl <- ring
    lbl$.ev_label_r <- lbl$arc_r1_var + label_gap
    lbl$.ev_hjust <- (1 - sin(lbl$.ev_mid_rad)) / 2
    lbl$.ev_vjust <- (1 - cos(lbl$.ev_mid_rad)) / 2
    lbl$.ev_lbl_x <- lbl$.ev_label_r * sin(lbl$.ev_mid_rad)
    lbl$.ev_lbl_y <- lbl$.ev_label_r * cos(lbl$.ev_mid_rad)
    lbl$.ev_lead_x <- (lbl$arc_r1_var + 0.1) * sin(lbl$.ev_mid_rad)
    lbl$.ev_lead_y <- (lbl$arc_r1_var + 0.1) * cos(lbl$.ev_mid_rad)
    lbl$.ev_lead_ex <- (lbl$.ev_label_r - 0.05) * sin(lbl$.ev_mid_rad)
    lbl$.ev_lead_ey <- (lbl$.ev_label_r - 0.05) * cos(lbl$.ev_mid_rad)
    lbl$.ev_lab_fill <- ifelse(lbl$.ev_is_up, pal$up, pal$down)
    max_r <- max(lbl$.ev_label_r) + label_headroom
    # Label boxes grow outward from their anchor; estimate their size (about
    # 0.05 units per character and 0.1 per line for each point of label_size)
    # so side labels are not clipped.
    lines <- strsplit(lbl$.ev_clean_label, "\n", fixed = TRUE)
    box_w <- vapply(lines, function(l) max(nchar(l)), numeric(1)) * label_size * 0.05
    box_h <- lengths(lines) * label_size * 0.1
    x_half <- max(
      max_r, abs(lbl$.ev_lbl_x - lbl$.ev_hjust * box_w) + label_headroom,
      abs(lbl$.ev_lbl_x + (1 - lbl$.ev_hjust) * box_w) + label_headroom
    )
    y_half <- max(
      max_r, abs(lbl$.ev_lbl_y - lbl$.ev_vjust * box_h) + label_headroom,
      abs(lbl$.ev_lbl_y + (1 - lbl$.ev_vjust) * box_h) + label_headroom
    )

    p <- p +
      ggplot2::geom_segment(
        data = lbl,
        ggplot2::aes(
          x = .data$.ev_lead_x, y = .data$.ev_lead_y,
          xend = .data$.ev_lead_ex, yend = .data$.ev_lead_ey
        ),
        linewidth = 0.4, colour = "grey35", inherit.aes = FALSE
      ) +
      ggplot2::geom_label(
        data = lbl,
        ggplot2::aes(
          x = .data$.ev_lbl_x, y = .data$.ev_lbl_y,
          label = .data$.ev_clean_label, fill = NULL
        ),
        hjust = lbl$.ev_hjust, vjust = lbl$.ev_vjust,
        fill = lbl$.ev_lab_fill, colour = "white", fontface = "bold",
        size = label_size, label.padding = ggplot2::unit(2, "pt"),
        label.r = ggplot2::unit(1.5, "pt"), lineheight = 0.85,
        inherit.aes = FALSE
      )
  }

  if (nrow(lab_pts) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = lab_pts,
      ggplot2::aes(.data$x_plot, .data$y_plot, label = .data$label_text),
      size = label_size, max.overlaps = 50, force = 4, seed = 42,
      inherit.aes = FALSE
    )
  }

  p +
    ggplot2::labs(
      title = title %||% "",
      subtitle = subtitle %||% NULL,
      tag = tag %||% NULL
    ) +
    ggplot2::coord_fixed(
      xlim = c(-1, 1) * ((x_half %||% max_r) + 0.1),
      ylim = c(-1, 1) * ((y_half %||% max_r) + 0.1), clip = "off"
    ) +
    ggplot2::theme_void(
      base_size = theme$base_size, base_family = theme$base_family
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 0)
      ),
      plot.subtitle = ggplot2::element_text(
        face = "bold.italic", hjust = 0.5,
        colour = "grey30",
        size = ggplot2::rel(0.65),
        margin = ggplot2::margin(t = 0, b = 0)
      ),
      plot.tag = ggplot2::element_text(face = "bold"),
      plot.tag.position = c(0.02, 0.99),
      plot.margin = ggplot2::margin(1, 1, 1, 1, "mm"),
      legend.position = "right",
      legend.key.height = ggplot2::unit(8, "mm"),
      legend.key.width = ggplot2::unit(2, "mm")
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(direction = "vertical"))
}
