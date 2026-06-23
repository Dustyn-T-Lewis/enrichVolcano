# Column resolution.
#
# Public functions take user-friendly column-name arguments and the surviving
# helpers downstream of `resolve_*()` use a list of resolved names. Keeping
# the lookup in one place avoids per-call grumbling about "what column was
# that again?".

#' Resolve volcano-frame columns
#'
#' @param volc_df Tidy DA tibble.
#' @param gene_col,logfc_col,pval_col,padj_col Column names.
#' @return Named list with the requested columns and a `has_padj` flag.
#' @keywords internal
#' @noRd
resolve_volc_cols <- function(volc_df, gene_col, logfc_col, pval_col, padj_col,
                              df_name = "volc_df") {
  for (cn in c(gene_col, logfc_col, pval_col)) {
    if (!cn %in% names(volc_df)) {
      ev_abort_missing_column(volc_df, cn, cn, df_name)
    }
  }
  list(
    gene     = gene_col,
    logfc    = logfc_col,
    pval     = pval_col,
    padj     = padj_col,
    has_padj = padj_col %in% names(volc_df)
  )
}

#' Resolve enrichment-frame columns
#'
#' @param enrich_df Tidy enrichment tibble.
#' @param term_col,nes_col,padj_col,size_col Column names.
#' @param magnitude One of `"neg_log_padj"`, `"size"`.
#' @return Named list with the requested columns.
#' @keywords internal
#' @noRd
resolve_enrich_cols <- function(enrich_df, term_col, nes_col, padj_col,
                                size_col, magnitude,
                                df_name = "enrich_df") {
  for (cn in c(term_col, padj_col, nes_col)) {
    if (!cn %in% names(enrich_df)) {
      ev_abort_missing_column(enrich_df, cn, cn, df_name)
    }
  }
  if (identical(magnitude, "size") && !size_col %in% names(enrich_df)) {
    ev_abort_missing_column(enrich_df, size_col, "size_col", df_name)
  }
  list(
    term = term_col,
    nes  = nes_col,
    padj = padj_col,
    size = size_col
  )
}

#' Resolve the leading-edge gene column (Q4 auto-detect)
#'
#' Walks `leading_edge` (`;`-string), `leadingEdge` (list-col),
#' `core_enrichment` (`/`-string), `Genes` (`;`-string). First match wins.
#' Emits one `cli_alert_info` naming the column + parser. No match -> `NULL`
#' (caller silently turns tick lines off).
#'
#' @param enrich_df Tidy enrichment tibble.
#' @param genes_col Optional override; if supplied, must exist in `enrich_df`.
#' @param genes_sep Optional separator when the column is a string.
#' @return A list of character vectors (one per row) or `NULL`.
#' @keywords internal
#' @noRd
resolve_genes_col <- function(enrich_df, genes_col, genes_sep,
                              df_name = "enrich_df") {
  if (!is.null(genes_col)) {
    if (!genes_col %in% names(enrich_df)) {
      ev_abort_missing_column(enrich_df, genes_col, "genes_col", df_name)
    }
    col <- enrich_df[[genes_col]]
    if (is.list(col)) {
      return(col)
    }
    return(strsplit(as.character(col), genes_sep %||% ";", fixed = TRUE))
  }
  candidates <- list(
    list(name = "leading_edge", sep = ";"),
    list(name = "leadingEdge", sep = NULL),
    list(name = "core_enrichment", sep = "/"),
    list(name = "Genes", sep = ";")
  )
  for (c in candidates) {
    if (c$name %in% names(enrich_df)) {
      col <- enrich_df[[c$name]]
      ev_inform("Reading tick lines from {.field {c$name}}.",
        class = "enrichVolcano_tick_resolution"
      )
      if (is.list(col)) {
        return(col)
      }
      return(strsplit(as.character(col), c$sep, fixed = TRUE))
    }
  }
  ev_inform("No tick-line column found; tick lines off.",
    class = "enrichVolcano_tick_resolution"
  )
  NULL
}
