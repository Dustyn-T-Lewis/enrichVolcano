# Volcano-in-ring composite: a single Cartesian ggplot with the volcano at the
# centre and enrichment pathways wrapped around it as an outer ring. Ported from
# the A_Proteomics_Analysis shared volcano_ring.R design. No coord_polar(); arcs
# are drawn with ggforce::geom_arc_bar() in a coord_fixed() panel.

# ---- label cleaning ----

#' Clean a pathway name for ring display
#'
#' Strips the database prefix, expands common acronyms, title-cases, wraps,
#' and applies MitoCarta-hierarchy shortening for `MITOCARTA_` pathways.
#'
#' @param name Character vector of raw pathway names.
#' @return Character vector of cleaned, wrapped labels.
#' @keywords internal
#' @noRd
ev_clean_label <- function(name) {
  if (length(name) > 1) {
    return(vapply(name, ev_clean_label, character(1), USE.NAMES = FALSE))
  }
  if (is.na(name) || !nzchar(name)) return(name)

  if (startsWith(name, "MITOCARTA_")) {
    return(ev_clean_label_mitocarta(name))
  }

  out <- name |>
    sub("^HALLMARK_", "", x = _) |>
    sub("^REACTOME_", "", x = _) |>
    sub("^GOSLIM_", "", x = _) |>
    sub("^KEGG_(MEDICUS|LEGACY)_", "", x = _) |>
    sub("^GO(BP|CC|MF)_", "", x = _) |>
    sub("^WP_", "", x = _) |>
    sub("^CORUM_", "", x = _) |>
    gsub("_", " ", x = _) |>
    tolower() |>
    tools::toTitleCase()

  out <- ev_expand_acronyms(out)
  out <- ev_shorten_phrases(out)
  out <- stringr::str_wrap(out, width = 15)
  trimws(out)
}

#' @keywords internal
#' @noRd
ev_expand_acronyms <- function(x) {
  reps <- c(
    "\\bTca\\b" = "TCA", "\\bMapk(\\d?)\\b" = "MAPK\\1",
    "\\bPtk(\\d?)\\b" = "PTK\\1", "\\bRhoa\\b" = "RhoA",
    "\\bRhob\\b" = "RhoB", "\\bRhoc\\b" = "RhoC",
    "\\bGtpase\\b" = "GTPase", "\\bGtp\\b" = "GTP", "\\bGdp\\b" = "GDP",
    "\\bMrna\\b" = "mRNA", "\\bRna\\b" = "RNA", "\\bDna\\b" = "DNA",
    "\\bAtp\\b" = "ATP", "\\bNadh\\b" = "NADH", "\\bFadh\\b" = "FADH",
    "\\bRos\\b" = "ROS", "\\bIfn\\b" = "IFN", "\\bIl-?(\\d+)\\b" = "IL\\1",
    "\\bPi3k\\b" = "PI3K", "\\bAkt\\b" = "AKT", "\\bMtor\\b" = "mTOR",
    "\\bMtorc1\\b" = "MTORC1", "\\bE2f\\b" = "E2F", "\\bMyc\\b" = "MYC",
    "\\bKras\\b" = "KRAS", "\\bTnfa\\b" = "TNFA", "\\bNfkb\\b" = "NFkB",
    "\\bG2m\\b" = "G2M", "\\bUv\\b" = "UV", "\\bWnt\\b" = "WNT",
    "\\bStat(\\d)\\b" = "STAT\\1", "\\bJak\\b" = "JAK", "\\bTgf\\b" = "TGF",
    "Trna " = "tRNA "
  )
  for (pat in names(reps)) x <- gsub(pat, reps[[pat]], x, perl = TRUE)
  x
}

#' @keywords internal
#' @noRd
ev_shorten_phrases <- function(x) {
  reps <- c(
    "Unfolded Protein Response" = "UPR", "Fatty Acid" = "FA",
    "Amino Acid" = "AA", "Generation Of" = "Gen. of", " And " = " & ",
    "Mitochondrion" = "Mito.", "Mitochondrial" = "Mito.",
    "Ubiquinone" = "UQ", "Organization" = "Org.",
    "Cytoskeleton" = "Cytoskel.", "Microtubule" = "MT",
    "Respiratory" = "Resp.", "Electron Transport" = "ETC",
    "Ubiquitin Dependent" = "Ub-Dep.", "Phosphorylation" = "Phosph.",
    "Modification" = "Mod.", "Intracellular" = "Intracell.",
    "Regulation Of" = "Reg.", "Signaling Pathway" = "Signaling",
    "Biosynthetic Process" = "Biosynthesis",
    "Catabolic Process" = "Catabolism", "Metabolic Process" = "Metabolism",
    "Response To" = "Resp. to", "Extracellular Matrix" = "ECM",
    "Epithelial Mesenchymal Transition" = "EMT"
  )
  for (pat in names(reps)) x <- gsub(pat, reps[[pat]], x, fixed = TRUE)
  x
}

#' @keywords internal
#' @noRd
ev_clean_label_mitocarta <- function(name) {
  n <- sub("^MITOCARTA_", "", name)
  parts <- strsplit(n, "__|>", perl = TRUE)[[1]]
  if (length(parts) == 1) parts <- strsplit(n, "_")[[1]]
  leaf <- utils::tail(parts, 1)
  leaf <- gsub("_", " ", leaf)
  leaf <- tools::toTitleCase(tolower(leaf))
  reps <- c(
    "\\bOxphos\\b" = "OXPHOS", "\\bTca\\b" = "TCA", "\\bAa\\b" = "AA",
    "\\bFa\\b" = "FA", "\\bRos\\b" = "ROS", "\\bImm\\b" = "IMM",
    "\\bOmm\\b" = "OMM", "\\bIms\\b" = "IMS",
    "\\bCi Subunits\\b" = "Complex I", "\\bCii Subunits\\b" = "Complex II",
    "\\bCiii Subunits\\b" = "Complex III", "\\bCiv Subunits\\b" = "Complex IV",
    "\\bCv Subunits\\b" = "Complex V", "\\bMitochondrial\\b" = "Mito.",
    "\\bMitochondrion\\b" = "Mito.", "\\bAmino Acid\\b" = "AA",
    "\\bFatty Acid\\b" = "FA"
  )
  for (pat in names(reps)) leaf <- gsub(pat, reps[[pat]], leaf, perl = TRUE)
  leaf <- trimws(gsub("\\s+", " ", leaf))
  stringr::str_wrap(leaf, width = 15)
}

# ---- ring geometry ----

# NES fill gradient (blue -> white -> red), shared by arcs and legend.
ev_nes_colours <- c("#08306B", "#4393C3", "white", "#D6604D", "#67000D")
ev_nes_values  <- c(-3, -1.5, 0, 1.5, 3)

#' Compute arc geometry for the enrichment ring
#'
#' Up-regulated pathways cluster at 90 deg (right), down at 270 deg (left);
#' a single term gets a compact 30 deg arc; arc outer radius encodes
#' `-log10(padj)`.
#' @keywords internal
#' @noRd
ev_ring_geometry <- function(enrich_df,
                             gap_intra = 3, gap_split = 8,
                             arc_r0 = 4.8, min_height = 0.05, max_height = 1.6) {
  if (nrow(enrich_df) == 0) return(enrich_df)

  ring <- enrich_df
  ring <- ring[!is.na(ring$pathway), , drop = FALSE]
  # fgsea sets padj/NES to NA for pathways with unbalanced stats; treat those
  # as non-significant (padj 1) with neutral NES so geometry stays finite.
  ring$padj[is.na(ring$padj)] <- 1
  ring$NES[is.na(ring$NES)] <- 0
  if (nrow(ring) == 0) return(ring)
  ring$clean_label <- ev_clean_label(ring$pathway)
  ring$gene_list <- strsplit(ring$leading_edge, ";", fixed = TRUE)
  ring$is_up <- !is.na(ring$NES) & ring$NES > 0

  up <- ring[ring$is_up, , drop = FALSE]
  dn <- ring[!ring$is_up, , drop = FALSE]
  up <- up[order(up$padj), , drop = FALSE]
  dn <- dn[order(dn$padj), , drop = FALSE]
  n_up <- nrow(up)
  n_dn <- nrow(dn)

  scale_h <- function(df) {
    if (nrow(df) == 0) return(numeric(0))
    if (nrow(df) == 1) return(arc_r0 + (min_height + max_height) / 2)
    neg_lp <- -log10(pmax(df$padj, .Machine$double.xmin))
    rng <- range(neg_lp)
    if (diff(rng) <= 0) {
      return(rep(arc_r0 + (min_height + max_height) / 2, nrow(df)))
    }
    scaled <- (neg_lp - rng[1]) / diff(rng)
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
  ring$mid_deg <- (ring$start_deg + ring$end_deg) / 2
  ring$start_rad <- ring$start_deg * pi / 180
  ring$end_rad <- ring$end_deg * pi / 180
  ring$mid_rad <- ring$mid_deg * pi / 180
  ring$term_idx <- seq_len(nrow(ring))
  ring
}

#' Per-gene radial ticks for each pathway arc
#' @keywords internal
#' @noRd
ev_tick_data <- function(ring_data, volc_df, tick_r0 = 4.4, tick_r1 = 4.8) {
  if (nrow(ring_data) == 0) return(tibble::tibble())
  gene_lfc <- volc_df[!is.na(volc_df$logFC), c("gene", "logFC"), drop = FALSE]
  gene_lfc <- gene_lfc[!duplicated(gene_lfc$gene), , drop = FALSE]
  pad <- 0.5 * pi / 180

  purrr::map_dfr(seq_len(nrow(ring_data)), function(i) {
    row <- ring_data[i, ]
    genes <- intersect(row$gene_list[[1]], gene_lfc$gene)
    if (length(genes) == 0) return(tibble::tibble())
    a0 <- row$start_rad + pad
    a1 <- row$end_rad - pad
    if (a1 <= a0) a1 <- a0 + pad
    ang <- seq(a0, a1, length.out = length(genes))
    lfc <- gene_lfc$logFC[match(genes, gene_lfc$gene)]
    tibble::tibble(
      gene = genes,
      logFC = lfc,
      direction = ifelse(lfc > 0, "up", "down"),
      x0 = tick_r0 * sin(ang), y0 = tick_r0 * cos(ang),
      x1 = tick_r1 * sin(ang), y1 = tick_r1 * cos(ang)
    )
  })
}

# ---- the composite ----

#' Build a volcano-in-ring composite for one contrast
#'
#' @param volc_df Tidy long tibble for a single contrast with `gene`, `logFC`,
#'   `P.Value`, and `pi_eq2` columns.
#' @param enrich_df Enrichment tibble for the same contrast (`pathway`, `padj`,
#'   `NES`, `size`, `leading_edge`).
#' @param title Plot title (defaults to the contrast name).
#' @param volcano_radius Inner volcano radius.
#' @param p_threshold Significance cutoff applied to `pi_eq2` for up/down calls.
#' @param point_size,point_alpha Volcano point aesthetics.
#' @param label_size Pathway-label text size.
#' @param disc_color Optional fill for a tinted central disc (default NULL = none).
#' @param theme Output of `ev_theme()`.
#' @return A ggplot object (class `c("enrichVolcano", ...)`).
#' @export
ev_volcano_ring <- function(volc_df, enrich_df, title = NULL,
                            volcano_radius = 3.5, p_threshold = 0.05,
                            disc_color = NULL,
                            point_size = 0.9, point_alpha = 0.6,
                            label_size = 2.6, theme = ev_theme()) {
  pal <- theme$palette
  vr <- volcano_radius * 0.92

  v <- volc_df[!is.na(volc_df$logFC) & !is.na(volc_df$P.Value), , drop = FALSE]
  v$neg_log10p <- -log10(v$P.Value)
  v <- v[is.finite(v$neg_log10p), , drop = FALSE]
  score <- if ("pi_eq2" %in% colnames(v)) v$pi_eq2 else v$P.Value
  v$direction <- ifelse(score < p_threshold & v$logFC > 0, "up",
                        ifelse(score < p_threshold & v$logFC < 0, "down", "ns"))
  n_up <- sum(v$direction == "up")
  n_down <- sum(v$direction == "down")

  x_max <- max(abs(v$logFC), na.rm = TRUE)
  y_max <- max(v$neg_log10p, na.rm = TRUE)
  if (!is.finite(x_max) || x_max == 0) x_max <- 1
  if (!is.finite(y_max) || y_max == 0) y_max <- 1
  v$x_plot <- v$logFC / x_max * vr
  v$y_plot <- v$neg_log10p / y_max * 2 * vr - vr
  v_ns <- v[v$direction == "ns", , drop = FALSE]
  v_sig <- v[v$direction != "ns", , drop = FALSE]

  ring <- ev_ring_geometry(enrich_df)
  ticks <- ev_tick_data(ring, v)

  p <- ggplot2::ggplot()

  if (!is.null(disc_color)) {
    disc <- data.frame(
      x = 4.4 * cos(seq(0, 2 * pi, length.out = 200)),
      y = 4.4 * sin(seq(0, 2 * pi, length.out = 200))
    )
    p <- p + ggplot2::geom_polygon(
      data = disc, ggplot2::aes(.data$x, .data$y),
      fill = disc_color, alpha = 0.12, colour = NA, inherit.aes = FALSE
    )
  }

  if (nrow(ring) > 0) {
    p <- p + ggforce::geom_arc_bar(
      data = ring,
      ggplot2::aes(x0 = 0, y0 = 0, r0 = 4.4, r = 4.8,
                   start = .data$start_rad, end = .data$end_rad),
      fill = "grey93", colour = "grey78", linewidth = 0.15,
      inherit.aes = FALSE
    )
  }

  p <- p +
    ggplot2::geom_point(
      data = v_ns, ggplot2::aes(.data$x_plot, .data$y_plot),
      colour = pal$ns, size = point_size * 0.7, alpha = point_alpha * 0.4
    ) +
    ggplot2::geom_point(
      data = v_sig,
      ggplot2::aes(.data$x_plot, .data$y_plot, colour = .data$direction),
      size = point_size, alpha = point_alpha
    ) +
    ggplot2::scale_colour_manual(
      values = c(up = pal$up, down = pal$down), guide = "none"
    ) +
    ggplot2::annotate("segment", x = 0, xend = 0, y = -vr, yend = vr * 0.96,
                      linewidth = 0.3, linetype = "dashed", colour = "grey50",
                      arrow = ggplot2::arrow(length = ggplot2::unit(1.2, "mm"),
                                             type = "closed")) +
    ggplot2::annotate("segment", x = -vr * 0.42, xend = vr * 0.42,
                      y = -vr, yend = -vr, linewidth = 0.3,
                      linetype = "dashed", colour = "grey50",
                      arrow = ggplot2::arrow(ends = "both",
                                             length = ggplot2::unit(1.2, "mm"),
                                             type = "closed")) +
    ggplot2::annotate("text", x = 0, y = vr * 1.05,
                      label = "-log10 p", size = 2.4, colour = "grey40") +
    ggplot2::annotate("text", x = vr * 0.45, y = -vr, label = "up",
                      size = 2.6, colour = pal$up,
                      fontface = "bold.italic", hjust = 0) +
    ggplot2::annotate("text", x = -vr * 0.45, y = -vr, label = "down",
                      size = 2.6, colour = pal$down,
                      fontface = "bold.italic", hjust = 1) +
    ggplot2::annotate("label", x = vr * 0.5, y = vr, label = n_up,
                      size = 3.2, colour = "white", fill = pal$up,
                      fontface = "bold", linewidth = 0.3) +
    ggplot2::annotate("label", x = -vr * 0.5, y = vr, label = n_down,
                      size = 3.2, colour = "white", fill = pal$down,
                      fontface = "bold", linewidth = 0.3)

  max_r <- 5.6
  if (nrow(ring) > 0) {
    if (nrow(ticks) > 0) {
      p <- p + ggplot2::geom_segment(
        data = ticks,
        ggplot2::aes(x = .data$x0, y = .data$y0, xend = .data$x1,
                     yend = .data$y1, colour = .data$direction),
        linewidth = 0.15, alpha = 0.7, inherit.aes = FALSE
      )
    }
    p <- p +
      ggforce::geom_arc_bar(
        data = ring,
        ggplot2::aes(x0 = 0, y0 = 0, r0 = 4.8, r = .data$arc_r1_var,
                     start = .data$start_rad, end = .data$end_rad,
                     fill = .data$NES),
        colour = "grey40", linewidth = 0.2, inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_gradientn(
        colours = ev_nes_colours,
        values = scales::rescale(ev_nes_values),
        limits = c(-3, 3), oob = scales::squish, name = "NES"
      )

    lbl <- ring
    lbl$label_r <- lbl$arc_r1_var + 1.4
    lbl$lbl_x <- lbl$label_r * sin(lbl$mid_rad)
    lbl$lbl_y <- lbl$label_r * cos(lbl$mid_rad)
    lbl$lead_x <- (lbl$arc_r1_var + 0.1) * sin(lbl$mid_rad)
    lbl$lead_y <- (lbl$arc_r1_var + 0.1) * cos(lbl$mid_rad)
    lbl$lead_ex <- (lbl$label_r - 0.3) * sin(lbl$mid_rad)
    lbl$lead_ey <- (lbl$label_r - 0.3) * cos(lbl$mid_rad)
    lbl$lab_fill <- ifelse(lbl$is_up, pal$up, pal$down)
    max_r <- max(lbl$label_r) + 1.2

    p <- p +
      ggplot2::geom_segment(
        data = lbl,
        ggplot2::aes(x = .data$lead_x, y = .data$lead_y,
                     xend = .data$lead_ex, yend = .data$lead_ey),
        linewidth = 0.4, colour = "grey35", inherit.aes = FALSE
      ) +
      ggplot2::geom_label(
        data = lbl,
        ggplot2::aes(x = .data$lbl_x, y = .data$lbl_y,
                     label = .data$clean_label, fill = NULL),
        fill = lbl$lab_fill, colour = "white", fontface = "bold",
        size = label_size, label.padding = ggplot2::unit(2, "pt"),
        label.r = ggplot2::unit(1.5, "pt"), lineheight = 0.85,
        inherit.aes = FALSE
      )
  }

  p <- p +
    ggplot2::labs(title = title %||% "") +
    ggplot2::coord_fixed(
      xlim = c(-(max_r + 0.6), max_r + 0.6),
      ylim = c(-(max_r + 0.6), max_r + 0.6), clip = "off"
    ) +
    ggplot2::theme_void(base_size = theme$base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.key.height = ggplot2::unit(10, "mm"),
      legend.key.width = ggplot2::unit(2.5, "mm")
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(direction = "vertical"))

  class(p) <- c("enrichVolcano", class(p))
  p
}
