#' Compute pi-score (Xiao 2014)
#'
#' @details
#' Two variants are supported, both from Xiao et al. 2014
#' *Bioinformatics* 30(6):801:
#' \itemize{
#'   \item \strong{Eq.2 (default):} \eqn{\pi = P^{|logFC|}}. Range \[0,1\],
#'         smaller = more significant; convention threshold 0.05.
#'   \item \strong{Eq.1 (alternative):} \eqn{\pi = |log_2 FC| \times -log_{10} P}.
#'         Unbounded, larger = more significant.
#' }
#'
#' @param data Tidy long tibble with `logFC` and `P.Value` columns.
#' @param variant `"eq2"` (default) or `"eq1"`.
#' @param p_col P-value column name.
#' @param logfc_col Log fold-change column name.
#' @param out_col Output column name (defaults to `pi_eq2` / `pi_eq1`).
#' @return Input tibble with one new column.
#' @references Xiao Y et al. (2014) *Bioinformatics* 30(6):801–807.
#' @export
pi_score <- function(data,
                        variant = c("eq2", "eq1"),
                        p_col = "P.Value",
                        logfc_col = "logFC",
                        out_col = NULL) {
  variant <- match.arg(variant)
  out_col <- out_col %||% paste0("pi_", variant)
  if (out_col %in% colnames(data)) return(data)
  p <- data[[p_col]]
  fc <- data[[logfc_col]]
  data[[out_col]] <- switch(variant,
    eq2 = p ^ abs(fc),
    eq1 = abs(fc) * -log10(p)
  )
  data
}
