#' Flag significance and direction
#'
#' Adds two columns to `data`: `sig` (logical) and `direction` (`"up"`,
#' `"down"`, or `"ns"`). The gate logic depends on which significance value
#' you choose, and the choice should match how that value is defined:
#' \itemize{
#'   \item `"adj.P.Val"` (default) or `"P.Value"`: threshold on the p-value
#'         AND a separate `|logFC|` gate. This is the reviewer-safe convention
#'         in differential-abundance proteomics.
#'   \item `"pi_eq2"`: threshold on pi-Eq.2 only. Pi-Eq.2 already folds in
#'         fold change (\eqn{\Pi = P^{|\log_2 \mathrm{FC}|}}), so a second
#'         `|logFC|` gate would double-count and exceed Xiao 2014's original
#'         proposal — the convention is to use Pi-Eq.2 directly as an
#'         adjusted-p-like value, NOT to BH-correct it (more stringent than
#'         nominal p, less conservative than BH).
#' }
#' Pi-Eq.1 is not accepted as a significance gate — it is unbounded with no
#' natural threshold; use it for ranking instead (`ev_enrich(rank_by = "signed_pi")`
#' or sort on `pi_eq1`).
#'
#' @param data Tidy long tibble with at minimum `logFC` and the chosen `value`
#'   column.
#' @param value One of `"adj.P.Val"`, `"P.Value"`, `"pi_eq2"`.
#' @param p_cut Significance threshold (default `0.05`).
#' @param fc_cut `|logFC|` gate, ignored when `value = "pi_eq2"`. Default
#'   `log2(1.5)` (~0.585).
#' @return Input tibble with `sig` (logical) and `direction` (character) columns
#'   appended or overwritten.
#' @references
#' Xiao Y et al. (2014) A novel significance score for gene identification from
#' RNA-seq data. \emph{Bioinformatics} 30(6):801-807.
#' @export
apply_gate <- function(data,
                       value = c("adj.P.Val", "P.Value", "pi_eq2"),
                       p_cut = 0.05,
                       fc_cut = log2(1.5)) {
  value <- match.arg(value)
  if (!"logFC" %in% colnames(data))
    ev_abort("{.var logFC} column required to assign direction.",
             class = "ev_missing_logfc")
  if (!value %in% colnames(data))
    ev_abort(c("Column {.var {value}} not found.",
               i = "Run the upstream step first (e.g. {.fn adjust_p} or {.fn pi_score})."),
             class = "ev_missing_gate_column")

  stat <- data[[value]]
  sig <- if (value == "pi_eq2") {
    stat < p_cut
  } else {
    stat < p_cut & abs(data$logFC) > fc_cut
  }
  sig[is.na(sig)] <- FALSE
  data$sig <- sig
  data$direction <- ifelse(!sig, "ns",
                           ifelse(data$logFC > 0, "up", "down"))
  data
}
