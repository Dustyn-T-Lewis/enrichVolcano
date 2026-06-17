#' Adjust p-values
#'
#' @details
#' \itemize{
#'   \item \strong{BH} (Benjamini-Hochberg, default): field-standard for proteomics.
#'   \item \strong{Bonferroni}: strict family-wise; use for confirmation.
#'   \item \strong{q-value} (Storey 2003): estimates pi_0, more power than BH;
#'         requires `qvalue` package.
#'   \item \strong{IHW} (Ignatiadis 2016): BH weighted by covariate;
#'         requires `IHW` package and a numeric `covariate` column.
#' }
#'
#' @param data Tidy long tibble.
#' @param method One of "BH", "bonferroni", "qvalue", "IHW".
#' @param covariate Column name (for IHW).
#' @param alpha Significance threshold (for IHW).
#' @param p_col Raw p-value column name.
#' @param by_contrast Logical; adjust within each contrast (default TRUE).
#' @return Tibble with `adj.P.Val` column (added or overwritten).
#' @export
adjust_p <- function(data,
                        method = c("BH", "bonferroni", "qvalue", "IHW"),
                        covariate = NULL,
                        alpha = 0.05,
                        p_col = "P.Value",
                        by_contrast = TRUE) {
  method <- match.arg(method)
  if (method == "IHW" && is.null(covariate)) {
    ev_abort("IHW requires a `covariate` column.",
             class = "ev_ihw_missing_covariate")
  }
  if (method == "qvalue" && !requireNamespace("qvalue", quietly = TRUE)) {
    ev_abort("qvalue method requires the {.pkg qvalue} package (Suggests).",
             class = "ev_missing_qvalue")
  }
  if (method == "IHW" && !requireNamespace("IHW", quietly = TRUE)) {
    ev_abort("IHW method requires the {.pkg IHW} package (Suggests).",
             class = "ev_missing_IHW")
  }

  adjust_one <- function(df) {
    p <- df[[p_col]]
    df$adj.P.Val <- switch(method,
      BH         = stats::p.adjust(p, method = "BH"),
      bonferroni = stats::p.adjust(p, method = "bonferroni"),
      qvalue     = tryCatch(
        qvalue::qvalue(p)$qvalues,
        error = function(e) {
          # qvalue's pi0 spline fails on few / narrow-range p-values (small
          # contrasts). Fall back to BH rather than surfacing a raw error.
          ev_warn(c("qvalue failed ({conditionMessage(e)}); falling back to BH.",
                    i = "Too few p-values or insufficient range to estimate pi0."),
                  class = "ev_qvalue_fallback")
          stats::p.adjust(p, method = "BH")
        }),
      IHW        = IHW::adj_pvalues(IHW::ihw(p, df[[covariate]], alpha = alpha))
    )
    df
  }

  if (by_contrast && "contrast" %in% colnames(data)) {
    data |>
      dplyr::group_by(.data$contrast) |>
      dplyr::group_modify(~ adjust_one(.x)) |>
      dplyr::ungroup()
  } else {
    adjust_one(data)
  }
}
