# Which rows of an enrichment a plot draws. Shared by volcano_ring() and
# nes_scatter() so both read the same view of one object.

check_enrichment <- function(x, arg = "enrichment") {
  if (!S7::S7_inherits(x, enrichment)) {
    ev_abort(
      c(
        "{.arg {arg}} must be an {.cls enrichment}, not {.cls {class(x)[1]}}.",
        i = "Build one with {.fn as_enrichment}."
      ),
      class = "enrichVolcano_input_error"
    )
  }
}

contrast_rows <- function(results, contrast) {
  available <- unique(results$contrast)
  if (is.null(contrast)) {
    if (length(available) > 1) {
      ev_abort(
        "This object holds {length(available)} contrasts; choose one with {.arg contrast}: {.val {available}}.",
        class = "enrichVolcano_input_error"
      )
    }
    contrast <- available
  }
  if (!isTRUE(contrast %in% available)) {
    ev_abort(
      "No contrast {.val {contrast}} in this object; it holds {.val {available}}.",
      class = "enrichVolcano_input_error"
    )
  }
  results[results$contrast == contrast, , drop = FALSE]
}

filter_view <- function(results, databases, collapse) {
  if (!is.null(databases)) {
    present <- unique(results$database[!is.na(results$database)])
    if (length(present) == 0) {
      ev_inform(
        "These results carry no database labels, so {.arg databases} was not applied.",
        class = "enrichVolcano_databases_skipped"
      )
    } else if (!any(databases %in% present)) {
      ev_abort(
        c("None of {.val {databases}} is in these results.", i = "Available: {.val {present}}."),
        class = "enrichVolcano_input_error"
      )
    } else {
      results <- results[results$database %in% databases, , drop = FALSE]
    }
  }
  if (collapse && "dedup_status" %in% names(results)) {
    results <- results[is.na(results$dedup_status) | results$dedup_status != "redundant", , drop = FALSE]
  }
  results
}

ring_terms <- function(results, term_threshold, n_terms, terms) {
  if (!is.null(terms)) {
    missing <- setdiff(terms, results$term)
    if (length(missing) > 0) {
      ev_abort("{.arg terms} not found in this contrast: {.val {missing}}.",
        class = "enrichVolcano_input_error"
      )
    }
    return(results[results$term %in% terms, , drop = FALSE])
  }
  sig <- results[!is.na(results$padj) & results$padj < term_threshold, , drop = FALSE]
  sig <- sig[order(sig$padj), , drop = FALSE]
  rank_in_direction <- stats::ave(seq_len(nrow(sig)), sig$direction, FUN = seq_along)
  sig[rank_in_direction <= n_terms, , drop = FALSE]
}

default_magnitude <- function(magnitude, score_type) {
  if (is.null(magnitude)) {
    return(if (identical(score_type, "NES")) "neg_log_padj" else "size")
  }
  rlang::arg_match(magnitude, c("neg_log_padj", "size"))
}

default_score_limits <- function(score_type, results, term_threshold) {
  if (identical(score_type, "NES")) {
    return(c(-3, 3))
  }
  sig <- results$score[!is.na(results$padj) & results$padj < term_threshold]
  m <- max(abs(if (length(sig) > 0) sig else results$score), 0, na.rm = TRUE)
  if (m == 0) m <- 1
  c(-m, m)
}
