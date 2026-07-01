# Volcano-in-ring composite: a single Cartesian ggplot with the volcano at the
# centre and enrichment pathways wrapped around it as an outer ring. Arcs are
# drawn with ggforce::geom_arc_bar() in a coord_fixed() panel.

# Ring geometry.

ev_ring_r_inner <- 4.4
ev_ring_r_outer <- 4.8

#' @keywords internal
#' @noRd
ev_ring_geometry <- function(enrich_df, term_col, padj_col, nes_col,
                             magnitude_col, genes_list, order_by = "padj",
                             gap_intra = 3, gap_split = 8,
                             arc_r0 = ev_ring_r_outer,
                             min_height = 0.05, max_height = 1.6) {
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
  ring$.ev_clean_label <- ev_clean_label(ring[[term_col]])
  ring$.ev_genes <- genes_list
  ring$.ev_is_up <- ring[[nes_col]] > 0
  ring$.ev_magnitude <- magnitude_col

  up <- ring[ring$.ev_is_up, , drop = FALSE]
  dn <- ring[!ring$.ev_is_up, , drop = FALSE]
  ord_key <- function(df) {
    switch(order_by,
      padj = df[[padj_col]],
      nes  = -abs(df[[nes_col]])
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

  if (n_up > 0 && n_dn > 0) {
    total_gap <- 2 * gap_split + (n_up - 1) * gap_intra + (n_dn - 1) * gap_intra
    arc_w <- (360 - total_gap) / (n_up + n_dn)
    up_span <- n_up * arc_w + (n_up - 1) * gap_intra
    up$start_deg <- (90 - up_span / 2) + (seq_len(n_up) - 1) * (arc_w + gap_intra)
    up$end_deg <- up$start_deg + arc_w
    dn_span <- n_dn * arc_w + (n_dn - 1) * gap_intra
    dn_offset <- (n_dn - seq_len(n_dn)) * (arc_w + gap_intra)
    dn$start_deg <- (270 - dn_span / 2) + dn_offset
    dn$end_deg <- dn$start_deg + arc_w
  } else {
    only <- if (n_up > 0) up else dn
    if (nrow(only) == 1) {
      center <- if (n_up > 0) 90 else 270
      only$start_deg <- center - 15
      only$end_deg <- only$start_deg + 30
    } else {
      arc_w <- (360 - nrow(only) * gap_intra) / nrow(only)
      only$start_deg <- (seq_len(nrow(only)) - 1) * (arc_w + gap_intra)
      only$end_deg <- only$start_deg + arc_w
    }
    if (n_up > 0) up <- only else dn <- only
  }

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
ev_tick_data <- function(ring_data, volc_df, gene_col, logfc_col,
                         tick_r0 = ev_ring_r_inner, tick_r1 = ev_ring_r_outer) {
  if (nrow(ring_data) == 0) {
    return(data.frame())
  }
  gene_lfc <- volc_df[!is.na(volc_df[[logfc_col]]),
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


# Composite assembly.

#' Volcano-in-ring composite for one contrast
#'
#' Draws a differential-abundance volcano embedded in a ring of enrichment
#' terms. NES-coloured arcs sit around the volcano; tick lines drop from the
#' arcs to each pathway's leading-edge genes inside the volcano.
#'
#' Column-naming arguments default to limma + fgsea conventions. The tick
#' column is auto-detected from `leading_edge` / `leadingEdge` /
#' `core_enrichment` / `Genes` unless `genes_col` is supplied.
#'
#' @param volc_df Tidy DA tibble for a single contrast.
#' @param enrich_df Tidy enrichment tibble for the same contrast.
#' @param gene_col,logfc_col,pval_col,padj_col Column names in `volc_df`.
#' @param volc_sig_col Optional column in `volc_df` used to call point
#'   significance (e.g. a pi-value), decoupled from the enrichment `padj_col`.
#'   `NULL` falls back to `padj_col` then `pval_col`. The y-axis stays
#'   `-log10(pval_col)` regardless.
#' @param term_col,nes_col,size_col Column names in `enrich_df`. `padj_col`
#'   is reused for the enrichment-side adjusted-p column.
#' @param genes_col Optional column name in `enrich_df` carrying leading-edge
#'   genes. `NULL` triggers auto-detect (Q4).
#' @param genes_sep Separator for the leading-edge string. `NULL` defers to
#'   the auto-detected default.
#' @param p_threshold Cutoff applied to `padj_col` in `volc_df`.
#' @param logfc_threshold Effect-size cutoff; a point is called up/down only
#'   when `abs(logFC) >= logfc_threshold` as well as significant.
#' @param title,subtitle,tag Plot text.
#' @param volcano_radius Inner volcano radius.
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
#' @param disc_color Optional fill for a tinted central disc.
#' @param nes_limits Length-2 numeric or `NULL`; defaults to `c(-3, 3)`.
#' @param magnitude `"neg_log_padj"` (default) or `"size"`; controls arc
#'   thickness encoding.
#' @param arc_order Angular order of arcs within each up/down half:
#'   `"padj"` (default, lowest FDR first) or `"nes"` (strongest `abs(NES)`
#'   first). The up/down split itself is always by NES sign.
#' @param arc_height_range Length-2 numeric `c(min, max)` for the shortest
#'   and tallest arc; widen it to exaggerate the magnitude encoding.
#' @param show_counts Draw the up/down significant-point count badges.
#' @param point_size,point_alpha Volcano point size (default 1.1) and opacity.
#' @param label_size Pathway-label text size.
#' @param label_gap Radial gap between each arc's outer top and its own label
#'   (default 0.6). Anchoring per-arc keeps every leader line the same short
#'   length regardless of arc height; the label box grows outward from this
#'   point so it never overlaps the arc. Widen it to lengthen all leaders.
#' @param count_size Text size of the up/down count badges (default 2.4).
#' @param count_x_mult,count_y_mult Badge position as a fraction of the volcano
#'   radius (default 0.7).
#' @param axis_size Text size of the `up`/`down`/`log2 FC`/`-log10 p` axis
#'   annotations (default 2.2).
#' @param label_mode,label_n,label_rank_by,label_genes Volcano point labels.
#' @param theme Output of `volcano_ring_theme()`.
#' @return A ggplot.
#' @export
volcano_ring <- function(volc_df, enrich_df,
                         gene_col = "gene",
                         logfc_col = "logFC",
                         pval_col = "P.Value",
                         padj_col = "padj",
                         volc_sig_col = NULL,
                         term_col = "pathway",
                         nes_col = "NES",
                         size_col = "size",
                         genes_col = NULL,
                         genes_sep = NULL,
                         p_threshold = 0.05,
                         logfc_threshold = 0,
                         title = NULL,
                         subtitle = NULL,
                         tag = NULL,
                         volcano_radius = 4.0,
                         ring_radius = 4.8,
                         ring_thickness = 0.55,
                         tick_width = 0.3,
                         label_headroom = 0.5,
                         disc_color = NULL,
                         nes_limits = NULL,
                         magnitude = c("neg_log_padj", "size"),
                         arc_order = c("padj", "nes"),
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
                         theme = volcano_ring_theme()) {
  magnitude <- match.arg(magnitude)
  arc_order <- match.arg(arc_order)
  label_mode <- match.arg(label_mode)
  label_rank_by <- match.arg(label_rank_by)
  ev_assert_colour(disc_color)
  validate_ring_geometry(ring_radius, volcano_radius, arc_height_range, label_headroom)

  if (!is.data.frame(volc_df)) {
    ev_abort(
      "{.arg volc_df} must be a data.frame, not {.cls {class(volc_df)[1]}}.",
      class = "enrichVolcano_input_error"
    )
  }
  if (!is.data.frame(enrich_df)) {
    ev_abort(
      "{.arg enrich_df} must be a data.frame, not {.cls {class(enrich_df)[1]}}.",
      class = "enrichVolcano_input_error"
    )
  }

  vcols <- resolve_volc_cols(volc_df, gene_col, logfc_col, pval_col, padj_col)
  if (!is.null(volc_sig_col) && !volc_sig_col %in% names(volc_df)) {
    ev_abort_missing_column(volc_df, volc_sig_col, "volc_sig_col", "volc_df")
  }
  ecols <- resolve_enrich_cols(
    enrich_df, term_col, nes_col, padj_col,
    size_col, magnitude
  )
  validate_volc_df(volc_df, vcols)
  validate_enrich_df(enrich_df, ecols)
  enrich_df <- dedup_by_term(enrich_df, ecols)
  enrich_df <- enrich_df[!is.na(enrich_df[[ecols$term]]), , drop = FALSE]
  has_padj <- vcols$has_padj

  pal <- theme$palette
  nes_limits <- nes_limits %||% theme$nes_limits %||% c(-3, 3)
  vr <- volcano_radius * 0.92
  ring_r0 <- ring_radius
  ring_r1 <- ring_radius + ring_thickness

  v <- volc_df[!is.na(volc_df[[logfc_col]]) & !is.na(volc_df[[pval_col]]), ,
    drop = FALSE
  ]
  v$.ev_neg_log10p <- -log10(v[[pval_col]])
  v <- v[is.finite(v$.ev_neg_log10p), , drop = FALSE]

  sig_score <- if (!is.null(volc_sig_col)) {
    v[[volc_sig_col]]
  } else if (has_padj) {
    v[[padj_col]]
  } else {
    v[[pval_col]]
  }
  sig <- sig_score < p_threshold
  lfc <- v[[logfc_col]]
  v$.ev_direction <- ifelse(sig & lfc >= logfc_threshold, "up",
    ifelse(sig & lfc <= -logfc_threshold, "down", "ns")
  )
  n_up <- sum(v$.ev_direction == "up")
  n_down <- sum(v$.ev_direction == "down")

  x_max <- max(abs(lfc), na.rm = TRUE)
  y_max <- max(v$.ev_neg_log10p, na.rm = TRUE)
  if (!is.finite(x_max) || x_max == 0) x_max <- 1
  if (!is.finite(y_max) || y_max == 0) y_max <- 1
  v$.ev_x_plot <- lfc / x_max * vr
  v$.ev_y_plot <- v$.ev_neg_log10p / y_max * 2 * vr - vr
  v_ns <- v[v$.ev_direction == "ns", , drop = FALSE]
  v_sig <- v[v$.ev_direction != "ns", , drop = FALSE]

  v_labels <- v
  names(v_labels)[names(v_labels) == gene_col] <- "gene"
  names(v_labels)[names(v_labels) == logfc_col] <- "logFC"
  v_labels$x_plot <- v$.ev_x_plot
  v_labels$y_plot <- v$.ev_y_plot
  lab_mode_old <- switch(label_mode,
    by_significance = "top_total",
    by_genes        = "explicit",
    label_mode
  )
  lab_pts <- ev_select_labels(
    v_labels,
    mode = lab_mode_old, n = label_n, rank_by = label_rank_by,
    genes = label_genes, p_col = pval_col,
    p_threshold = p_threshold, logfc_threshold = logfc_threshold
  )

  mag_vec <- switch(magnitude,
    neg_log_padj = -log10(pmax(enrich_df[[ecols$padj]], .Machine$double.xmin)),
    size         = as.numeric(enrich_df[[ecols$size]])
  )
  genes_list <- resolve_genes_col(enrich_df, genes_col, genes_sep)
  if (is.null(genes_list)) {
    genes_list <- replicate(nrow(enrich_df), character(0), simplify = FALSE)
  }
  ring <- ev_ring_geometry(enrich_df,
    term_col = term_col, padj_col = padj_col,
    nes_col = nes_col, magnitude_col = mag_vec,
    genes_list = genes_list, order_by = arc_order, arc_r0 = ring_r1,
    min_height = arc_height_range[1], max_height = arc_height_range[2]
  )
  ticks <- ev_tick_data(ring, v,
    gene_col = gene_col, logfc_col = logfc_col,
    tick_r0 = ring_r0, tick_r1 = ring_r1
  )

  p <- ggplot2::ggplot()

  if (!is.null(disc_color)) {
    disc <- data.frame(
      x = ring_r0 * cos(seq(0, 2 * pi, length.out = 200)),
      y = ring_r0 * sin(seq(0, 2 * pi, length.out = 200))
    )
    p <- p + ggplot2::geom_polygon(
      data = disc, ggplot2::aes(.data$x, .data$y),
      fill = disc_color, alpha = 0.12, colour = NA, inherit.aes = FALSE
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
          fill = .data[[nes_col]]
        ),
        colour = "grey40", linewidth = 0.2, inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_gradientn(
        colours = pal$nes_scale,
        values = scales::rescale(pal$nes_values),
        limits = nes_limits, oob = scales::squish, name = "NES"
      )

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
      xlim = c(-(max_r + 0.1), max_r + 0.1),
      ylim = c(-(max_r + 0.1), max_r + 0.1), clip = "off"
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
