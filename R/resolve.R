# Column resolution. Public functions take user-facing column-name arguments;
# helpers downstream of `resolve_*()` work from the resolved-name list these
# return, so the name lookup lives in one place.

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
