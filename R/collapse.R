ev_subset_preserve_attr <- function(df, idx,
                                    attrs = c("ev_pathways", "ev_stats", "ev_filter")) {
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
#'     \code{fgsea::collapsePathways} first on significant pathways
#'     (\code{padj < sig_threshold}), then Jaccard on the survivors.
#'     Matches the source pipeline that generated the YvO 2025 figures.
#'   \item \strong{jaccard_then_collapse}: Jaccard first, then
#'     collapsePathways. Equivalent to the deprecated \code{"both"}.
#'   \item \strong{jaccard}: \eqn{J(A,B) = |A \cap B| / |A \cup B|}.
#'     Above cutoff, keep representative (by \code{keep_by}). Data-independent.
#'   \item \strong{collapse_fgsea}: \code{fgsea::collapsePathways} only.
#'     Data-dependent; preferred for single-database fgsea output.
#'     Requires the \code{ev_pathways} attribute attached by
#'     \code{ev_enrich()}; if absent, warns and keeps all pathways for
#'     that group.
#'   \item \strong{both}: deprecated alias for \code{jaccard_then_collapse}.
#' }
#'
#' Non-significant rows (\code{padj >= sig_threshold} or NA) always retain
#' \code{dedup_kept = TRUE}. Set \code{sig_threshold = NA} to disable the
#' gate and dedup every row.
#'
#' @param enrich_result Tibble from `ev_enrich()`.
#' @param method One of \code{"collapse_then_jaccard"} (default),
#'   \code{"jaccard_then_collapse"}, \code{"jaccard"}, or
#'   \code{"collapse_fgsea"}. \code{"both"} is a deprecated alias for
#'   \code{"jaccard_then_collapse"}.
#' @param cutoff Similarity cutoff. Length 1 (shared by all `scope` stages) or
#'   one value per stage.
#' @param scope `"within_db"` (collapse inside each database), `"cross_db"`
#'   (collapse across databases of a contrast), or `"global"`. May be a vector
#'   of stages run in order, e.g. `c("within_db", "cross_db")`: each stage
#'   collapses only the survivors of the previous, so you can clean each
#'   database's own hierarchy first and then merge the same concept across
#'   databases (optionally with a different `cutoff` per stage). Note `cross_db`
#'   alone already compares every pathway pairwise within a contrast, so it
#'   subsumes `within_db`; staging mainly matters when you want a tighter
#'   within-database cutoff than the cross-database one.
#' @param sig_threshold Numeric in \[0, 1\] or NA. Only rows with
#'   \code{padj < sig_threshold} are considered for dedup; non-significant
#'   rows always retain \code{dedup_kept = TRUE}. NA disables the gate.
#'   Default 0.05.
#' @param collapse_pval_threshold `pval.threshold` forwarded to
#'   `fgsea::collapsePathways`. Used by every method that invokes the fgsea
#'   collapse step: `collapse_fgsea`, `jaccard_then_collapse`, and the
#'   default `collapse_then_jaccard`. Default 0.05.
#' @param keep_by Tie-breaking column ("padj", "size", "NES").
#' @param similarity Gene-set similarity metric for the Jaccard-family steps,
#'   computed on leading-edge genes (Merico 2010, EnrichmentMap):
#'   \code{"jaccard"} (default, \eqn{|A\cap B|/|A\cup B|}; EM cutoff 0.25),
#'   \code{"overlap"} (\eqn{|A\cap B|/\min(|A|,|B|)}, catches a small set
#'   contained in a larger one; EM cutoff 0.5), or \code{"combined"}
#'   (\eqn{w\cdot Jaccard + (1-w)\cdot Overlap}; the Cytoscape EnrichmentMap
#'   default, cutoff 0.375). \code{cutoff} is interpreted on the chosen metric.
#' @param combined_weight Weight \eqn{w} on the Jaccard term for
#'   \code{similarity = "combined"} (default 0.5, per EnrichmentMap). The
#'   overlap term gets \eqn{1 - w}.
#' @return Tibble with `dedup_kept` logical column (TRUE = retained).
#' @export
ev_collapse <- function(enrich_result,
                        method = c("collapse_then_jaccard",
                                   "jaccard_then_collapse",
                                   "jaccard",
                                   "collapse_fgsea",
                                   "both"),
                        cutoff = 0.5,
                        scope = "within_db",
                        sig_threshold = 0.05,
                        collapse_pval_threshold = 0.05,
                        keep_by = c("padj", "size", "NES"),
                        similarity = c("jaccard", "overlap", "combined"),
                        combined_weight = 0.5) {
  method <- match.arg(method)
  keep_by <- match.arg(keep_by)
  similarity <- match.arg(similarity)
  # `scope` may be a vector of stages run in order (e.g. c("within_db",
  # "cross_db")): each stage collapses the survivors of the previous one.
  # `cutoff` is length 1 (shared) or one per stage.
  valid_scope <- c("within_db", "cross_db", "global")
  if (!all(scope %in% valid_scope)) {
    ev_abort("`scope` must be one or more of {.val {valid_scope}}.",
             class = "ev_bad_scope")
  }
  if (length(cutoff) == 1L) cutoff <- rep(cutoff, length(scope))
  if (length(cutoff) != length(scope)) {
    ev_abort("`cutoff` must be length 1 or match the number of `scope` stages.",
             class = "ev_bad_cutoff")
  }
  if (combined_weight < 0 || combined_weight > 1) {
    ev_abort("`combined_weight` must be in [0, 1]; got {combined_weight}.",
             class = "ev_bad_combined_weight")
  }
  if (identical(method, "both")) {
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = 'ev_collapse(method = "both")',
      with = 'ev_collapse(method = "jaccard_then_collapse")'
    )
    method <- "jaccard_then_collapse"
  }
  enrich_result$dedup_kept <- TRUE
  for (s in seq_along(scope)) {
    enrich_result <- ev_collapse_stage(
      enrich_result, method, cutoff[s], scope[s], sig_threshold,
      collapse_pval_threshold, keep_by, similarity, combined_weight)
  }
  enrich_result
}

# One dedup stage at a single scope. Respects any existing `dedup_kept` so
# stages chain: later scopes only consider rows still kept by earlier ones.
ev_collapse_stage <- function(enrich_result, method, cutoff, scope,
                              sig_threshold, collapse_pval_threshold,
                              keep_by, similarity, combined_weight) {
  group_cols <- switch(scope,
    within_db = c("contrast", "database"),
    cross_db  = "contrast",
    global    = character(0)
  )
  if (length(group_cols) > 0) {
    keys <- do.call(paste, c(enrich_result[group_cols], sep = "::"))
  } else {
    keys <- rep("global", nrow(enrich_result))
  }
  for (k in unique(keys)) {
    idx <- which(keys == k)
    if (length(idx) < 2) next
    sig_idx <- ev_sig_idx(idx, enrich_result$padj, sig_threshold)
    sig_idx <- sig_idx[enrich_result$dedup_kept[sig_idx]]   # survivors of prior stages
    if (length(sig_idx) < 2) next
    if (method %in% c("jaccard", "jaccard_then_collapse")) {
      enrich_result$dedup_kept[sig_idx] <- ev_setsim_dedup(
        ev_subset_preserve_attr(enrich_result, sig_idx), cutoff, keep_by,
        similarity, combined_weight
      )
    }
    if (method %in% c("collapse_fgsea", "jaccard_then_collapse")) {
      kept_sig <- sig_idx[enrich_result$dedup_kept[sig_idx]]
      if (length(kept_sig) >= 2) {
        # Pass only the Jaccard survivors (not the full sig_idx) and write
        # back positionally — mirroring the collapse_then_jaccard branch
        # below. The previous code ANDed the collapse result over sig_idx,
        # which silently dropped both members of a redundant cluster when
        # collapsePathways and Jaccard's keep_by tie-breaker disagreed on
        # which pathway represents the cluster.
        survive <- ev_collapse_fgsea(
          ev_subset_preserve_attr(enrich_result, kept_sig),
          collapse_pval_threshold
        )
        enrich_result$dedup_kept[kept_sig] <- survive
      }
    }
    if (identical(method, "collapse_then_jaccard")) {
      survive_collapse <- ev_collapse_fgsea(
        ev_subset_preserve_attr(enrich_result, sig_idx),
        collapse_pval_threshold
      )
      kept_after_collapse <- sig_idx[survive_collapse]
      if (length(kept_after_collapse) >= 2) {
        survive_jaccard <- ev_setsim_dedup(
          ev_subset_preserve_attr(enrich_result, kept_after_collapse),
          cutoff, keep_by, similarity, combined_weight
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

# Set-similarity dedup on leading-edge gene sets (Merico 2010 PMID:21085593).
#   jaccard  = |A∩B| / |A∪B|              symmetric overlap
#   overlap  = |A∩B| / min(|A|, |B|)      catches asymmetric containment
#   combined = w*Jaccard + (1-w)*Overlap  EnrichmentMap default (w = 0.5)
# Iterate best-first by keep_by; a later set is dropped if its similarity to a
# retained set exceeds cutoff.
ev_setsim_dedup <- function(df, cutoff, keep_by,
                            similarity = "jaccard", combined_weight = 0.5) {
  sets <- strsplit(df$leading_edge, ";", fixed = TRUE)
  n <- length(sets)
  kept <- rep(TRUE, n)
  # Rank best-first: smallest padj, but LARGEST size and largest |NES| (a strong
  # down-regulated term must outrank a weak up-regulated one — signed NES would
  # invert that).
  rank_key <- switch(keep_by,
    padj = df$padj,
    size = -df$size,
    NES  = -abs(df$NES)
  )
  ord <- order(rank_key)
  # Only collapse pathways of the SAME direction — an up- and a down-regulated
  # term that happen to share leading-edge genes are distinct biology.
  dir <- if ("direction" %in% names(df)) df$direction else rep(NA_character_, n)
  for (i in ord) {
    if (!kept[i]) next
    for (j in ord) {
      if (i == j || !kept[j]) next
      if (!is.na(dir[i]) && !is.na(dir[j]) && dir[i] != dir[j]) next
      inter <- length(intersect(sets[[i]], sets[[j]]))
      uni <- length(union(sets[[i]], sets[[j]]))
      min_n <- min(length(sets[[i]]), length(sets[[j]]))
      jac <- if (uni > 0) inter / uni else 0
      ovl <- if (min_n > 0) inter / min_n else 0
      sim <- switch(similarity,
        jaccard  = jac,
        overlap  = ovl,
        combined = combined_weight * jac + (1 - combined_weight) * ovl
      )
      if (sim > cutoff) kept[j] <- FALSE
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
  kept <- rep(TRUE, nrow(df))
  fg_mask <- df$mode == "fgsea"
  if (!any(fg_mask)) return(kept)

  # ev_enrich attaches ev_pathways keyed by database, each value itself a
  # named list of pathway -> gene character vectors. Test fixtures sometimes
  # attach a flat name-keyed list. Detect the shape and look up by
  # (database, pathway) when db-keyed so cross-DB pathway-name collisions
  # resolve to the correct gene set rather than to the first matching name.
  db_keyed <- length(pathway_list) > 0 &&
    all(vapply(pathway_list, is.list, logical(1)))
  resolve_pw <- function(db, pw) {
    if (db_keyed) pathway_list[[db]][[pw]] else pathway_list[[pw]]
  }

  # collapsePathways is RNG-sensitive and consumes a single per-contrast rank
  # vector. Partition by contrast so each subset is scored against the
  # matching stats; under scope = "global" with multiple contrasts the prior
  # implementation silently scored every pathway against contrast[1]'s ranks.
  fg_contrasts <- unique(df$contrast[fg_mask])
  for (ctr in fg_contrasts) {
    ctr_mask <- fg_mask & df$contrast == ctr
    ctr_rows <- df[ctr_mask, , drop = FALSE]
    stats <- stats_list[[ctr]]
    if (is.null(stats) || length(stats) == 0) {
      ev_warn(
        "No gene-level stats for contrast {.val {ctr}}; keeping all pathways.",
        class = "ev_collapse_no_pathways"
      )
      next
    }
    pw_for_fg <- Map(resolve_pw, ctr_rows$database, ctr_rows$pathway)
    has_paths <- !vapply(pw_for_fg, is.null, logical(1))
    if (!any(has_paths)) next
    eval_rows <- ctr_rows[has_paths, , drop = FALSE]
    pw_for_fg <- pw_for_fg[has_paths]
    names(pw_for_fg) <- eval_rows$pathway
    fg_input <- data.table::data.table(
      pathway     = eval_rows$pathway,
      pval        = eval_rows$pval,
      padj        = eval_rows$padj,
      log2err     = if ("log2err" %in% colnames(eval_rows)) eval_rows$log2err else NA_real_,
      ES          = if ("ES" %in% colnames(eval_rows)) eval_rows$ES else eval_rows$NES,
      NES         = eval_rows$NES,
      size        = eval_rows$size,
      leadingEdge = strsplit(eval_rows$leading_edge, ";", fixed = TRUE)
    )
    # collapsePathways re-runs fgsea on residual gene lists; pin the seed so
    # dedup is reproducible across calls.
    collapsed <- tryCatch(
      withr::with_seed(42, {
        fgsea::collapsePathways(
          fgseaRes       = fg_input,
          pathways       = pw_for_fg,
          stats          = stats,
          pval.threshold = pval_threshold
        )
      }),
      error = function(e) {
        ev_warn(
          "fgsea::collapsePathways failed: {conditionMessage(e)}; keeping all pathways.",
          class = "ev_collapse_fgsea_error"
        )
        NULL
      }
    )
    if (is.null(collapsed)) next
    eval_df_idx <- which(ctr_mask)[has_paths]
    kept[eval_df_idx] <- eval_rows$pathway %in% collapsed$mainPathways
  }
  kept
}
