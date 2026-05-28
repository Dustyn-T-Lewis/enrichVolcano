ev_subset_preserve_attr <- function(df, idx, attrs = c("ev_pathways", "ev_stats")) {
  out <- df[idx, , drop = FALSE]
  for (a in attrs) attr(out, a) <- attr(df, a)
  out
}

ev_sig_idx <- function(idx, padj, sig_threshold) {
  if (is.na(sig_threshold)) return(idx)
  idx[!is.na(padj[idx]) & padj[idx] < sig_threshold]
}

#' Deduplicate enrichment results
#'
#' @details
#' \itemize{
#'   \item \strong{collapse_then_jaccard} (default, since 0.2.0): runs
#'     `fgsea::collapsePathways` first on significant rows
#'     (`padj < sig_threshold`), then Jaccard on the survivors. Matches
#'     the source pipeline used for the YvO 2025 figures.
#'   \item \strong{jaccard_then_collapse}: Jaccard first, then
#'     collapsePathways. Equivalent to the deprecated `"both"`.
#'   \item \strong{jaccard}: \eqn{J(A,B) = |A \cap B| / |A \cup B|}.
#'     Above cutoff, keep representative (by `keep_by`). Data-independent.
#'   \item \strong{collapse_fgsea}: `fgsea::collapsePathways` — tests
#'     leading-edge dependency. Data-dependent; preferred for fgsea output.
#'     Requires the `ev_pathways` attribute attached by `ev_enrich()`; if
#'     absent, warns and keeps all pathways for that group.
#'   \item \strong{both}: deprecated alias for `jaccard_then_collapse`.
#' }
#'
#' @param enrich_result Tibble from `ev_enrich()`.
#' @param method One of `"collapse_then_jaccard"` (default), `"jaccard_then_collapse"`,
#'   `"jaccard"`, or `"collapse_fgsea"`. `"both"` is a deprecated alias for
#'   `"jaccard_then_collapse"`.
#' @param cutoff Jaccard similarity cutoff.
#' @param scope `"within_db"`, `"cross_db"`, or `"global"`.
#' @param sig_threshold Numeric in \[0, 1\] or NA. Only rows with
#'   `padj < sig_threshold` enter dedup. Non-significant rows always retain
#'   `dedup_kept = TRUE`. NA disables the gate. Default 0.05.
#' @param collapse_pval_threshold For `collapse_fgsea`.
#' @param keep_by Tie-breaking column ("padj", "size", "NES").
#' @return Tibble with `dedup_kept` logical column (TRUE = retained).
#' @export
ev_collapse <- function(enrich_result,
                        method = c("collapse_then_jaccard",
                                   "jaccard_then_collapse",
                                   "jaccard",
                                   "collapse_fgsea",
                                   "both"),
                        cutoff = 0.5,
                        scope = c("within_db", "cross_db", "global"),
                        sig_threshold = 0.05,
                        collapse_pval_threshold = 0.05,
                        keep_by = c("padj", "size", "NES")) {
  method <- match.arg(method)
  scope <- match.arg(scope)
  keep_by <- match.arg(keep_by)
  if (identical(method, "both")) {
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = 'ev_collapse(method = "both")',
      with = 'ev_collapse(method = "jaccard_then_collapse")'
    )
    method <- "jaccard_then_collapse"
  }
  group_cols <- switch(scope,
    within_db = c("contrast", "database"),
    cross_db  = "contrast",
    global    = character(0)
  )
  enrich_result$dedup_kept <- TRUE
  if (length(group_cols) > 0) {
    keys <- do.call(paste, c(enrich_result[group_cols], sep = "::"))
  } else {
    keys <- rep("global", nrow(enrich_result))
  }
  for (k in unique(keys)) {
    idx <- which(keys == k)
    if (length(idx) < 2) next
    sig_idx <- ev_sig_idx(idx, enrich_result$padj, sig_threshold)
    if (length(sig_idx) < 2) next
    if (method %in% c("jaccard", "jaccard_then_collapse")) {
      enrich_result$dedup_kept[sig_idx] <- ev_jaccard_dedup(
        ev_subset_preserve_attr(enrich_result, sig_idx), cutoff, keep_by
      )
    }
    if (method %in% c("collapse_fgsea", "jaccard_then_collapse")) {
      kept_sig <- sig_idx[enrich_result$dedup_kept[sig_idx]]
      if (length(kept_sig) >= 2) {
        enrich_result$dedup_kept[sig_idx] <- enrich_result$dedup_kept[sig_idx] &
          ev_collapse_fgsea(ev_subset_preserve_attr(enrich_result, sig_idx),
                            collapse_pval_threshold)
      }
    }
    if (identical(method, "collapse_then_jaccard")) {
      survive_collapse <- ev_collapse_fgsea(
        ev_subset_preserve_attr(enrich_result, sig_idx),
        collapse_pval_threshold
      )
      kept_after_collapse <- sig_idx[survive_collapse]
      if (length(kept_after_collapse) >= 2) {
        survive_jaccard <- ev_jaccard_dedup(
          ev_subset_preserve_attr(enrich_result, kept_after_collapse),
          cutoff, keep_by
        )
        survive_final <- rep(FALSE, length(sig_idx))
        survive_final[survive_collapse] <- survive_jaccard
        enrich_result$dedup_kept[sig_idx] <- survive_final
      } else {
        enrich_result$dedup_kept[sig_idx] <- survive_collapse
      }
    }
  }
  enrich_result
}

ev_jaccard_dedup <- function(df, cutoff, keep_by) {
  sets <- strsplit(df$leading_edge, ";", fixed = TRUE)
  n <- length(sets)
  kept <- rep(TRUE, n)
  ord <- order(df[[keep_by]], decreasing = keep_by != "padj")
  for (i in ord) {
    if (!kept[i]) next
    for (j in ord) {
      if (i == j || !kept[j]) next
      inter <- length(intersect(sets[[i]], sets[[j]]))
      uni <- length(union(sets[[i]], sets[[j]]))
      jac <- if (uni > 0) inter / uni else 0
      if (jac > cutoff) kept[j] <- FALSE
    }
  }
  kept
}

ev_collapse_fgsea <- function(df, pval_threshold) {
  pathway_list <- attr(df, "ev_pathways")
  stats_list   <- attr(df, "ev_stats")
  if (is.null(pathway_list) || is.null(stats_list)) {
    ev_warn(
      c("ev_collapse(method = 'collapse_fgsea') requires the 'ev_pathways' and",
        "'ev_stats' attributes from ev_enrich(); keeping all pathways for this group."),
      class = "ev_collapse_no_pathways"
    )
    return(rep(TRUE, nrow(df)))
  }
  fg_subset <- df[df$mode == "fgsea", , drop = FALSE]
  if (nrow(fg_subset) == 0) return(rep(TRUE, nrow(df)))
  ctr   <- fg_subset$contrast[1]
  stats <- stats_list[[ctr]]
  if (is.null(stats) || length(stats) == 0) {
    ev_warn(
      "No gene-level stats for contrast {.val {ctr}}; keeping all pathways.",
      class = "ev_collapse_no_pathways"
    )
    return(rep(TRUE, nrow(df)))
  }
  fg_input <- data.table::data.table(
    pathway     = fg_subset$pathway,
    pval        = fg_subset$pval,
    padj        = fg_subset$padj,
    log2err     = if ("log2err" %in% colnames(fg_subset)) fg_subset$log2err else NA_real_,
    ES          = if ("ES" %in% colnames(fg_subset)) fg_subset$ES else fg_subset$NES,
    NES         = fg_subset$NES,
    size        = fg_subset$size,
    leadingEdge = strsplit(fg_subset$leading_edge, ";", fixed = TRUE)
  )
  collapsed <- tryCatch(
    fgsea::collapsePathways(
      fgseaRes       = fg_input,
      pathways       = pathway_list[fg_subset$pathway],
      stats          = stats,
      pval.threshold = pval_threshold
    ),
    error = function(e) {
      ev_warn(
        "fgsea::collapsePathways failed: {conditionMessage(e)}; keeping all pathways.",
        class = "ev_collapse_fgsea_error"
      )
      NULL
    }
  )
  if (is.null(collapsed)) return(rep(TRUE, nrow(df)))
  main   <- collapsed$mainPathways
  kept   <- rep(TRUE, nrow(df))
  fg_idx <- which(df$mode == "fgsea")
  kept[fg_idx] <- df$pathway[fg_idx] %in% main
  kept
}
