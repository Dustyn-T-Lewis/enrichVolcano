# Shared rule-driven label selection for the standalone volcano and the
# composite. Default mode "none" labels nothing.

#' Label text for a row: symbol, else accession, else gene.
#' @keywords internal
#' @noRd
ev_label_text <- function(df) {
  symbol  <- if ("symbol" %in% names(df)) df$symbol else NA_character_
  uniprot <- if ("uniprot" %in% names(df)) df$uniprot else NA_character_
  gene    <- if ("gene" %in% names(df)) df$gene else NA_character_
  out <- ifelse(!is.na(symbol) & nzchar(symbol), symbol,
                ifelse(!is.na(uniprot) & nzchar(uniprot), uniprot, gene))
  as.character(out)
}

#' Select which rows to label on a volcano
#'
#' @param df Single-contrast tibble with `logFC` and the significance column.
#' @param mode One of `"none"`, `"top_per_direction"`, `"top_total"`,
#'   `"all_significant"`, `"explicit"`.
#' @param n Cap for the top-* modes (per direction for `top_per_direction`).
#' @param rank_by `"significance"` (smaller `p_col` first) or `"logfc"`
#'   (larger `|logFC|` first).
#' @param genes For `"explicit"`: matched against symbol, accession, or gene.
#' @param p_col Significance column (e.g. `"pi_eq2"`, `"P.Value"`,
#'   `"adj.P.Val"`).
#' @param p_threshold,logfc_threshold Cutoffs defining the candidate pool.
#' @return Sub-tibble of rows to label, with a `label_text` column.
#' @keywords internal
#' @noRd
ev_select_labels <- function(df, mode, n, rank_by, genes,
                             p_col, p_threshold, logfc_threshold) {
  df$label_text <- ev_label_text(df)
  if (mode == "explicit") {
    keep <- df$label_text %in% genes |
      ("symbol" %in% names(df) & df$symbol %in% genes) |
      ("uniprot" %in% names(df) & df$uniprot %in% genes) |
      df$gene %in% genes
    return(df[keep & !is.na(keep), , drop = FALSE])
  }
  if (mode == "none") return(df[0, , drop = FALSE])

  sig <- df[[p_col]] < p_threshold & abs(df$logFC) >= logfc_threshold
  pool <- df[sig & !is.na(sig), , drop = FALSE]
  if (nrow(pool) == 0) return(pool)

  ord <- if (rank_by == "logfc") {
    order(-abs(pool$logFC))
  } else {
    order(pool[[p_col]])
  }
  pool <- pool[ord, , drop = FALSE]

  switch(mode,
    all_significant = pool,
    top_total = utils::head(pool, n),
    top_per_direction = {
      up <- pool[pool$logFC > 0, , drop = FALSE]
      dn <- pool[pool$logFC < 0, , drop = FALSE]
      rbind(utils::head(up, n), utils::head(dn, n))
    },
    ev_abort("Unknown label mode {.val {mode}}.", class = "ev_bad_label_mode")
  )
}
