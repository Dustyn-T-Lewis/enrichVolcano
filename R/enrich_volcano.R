#' Build a composite volcano + enrichment ring plot
#'
#' Hero function for `enrichVolcano`. Runs the full pipeline:
#' validate input -> compute pi-score -> adjust p-values -> run enrichment ->
#' deduplicate -> build volcano + ring panels -> compose patchwork.
#'
#' @param data Differential abundance results. Accepts tidy long, wide-suffix,
#'   DEP, proteoDA, limma, DESeq2, edgeR, MSstats, proDA, DEqMS, MaxQuant,
#'   or Perseus formats - `ev_validate()` auto-detects.
#' @param contrast Character; one or many contrast names.
#' @param species Character; "human" (default), "mouse", "rat", "zebrafish",
#'   "fly", "yeast", "pig".
#' @param databases Character vector of registered DB names (see
#'   `list_databases()`), or a named list of pathway lists.
#' @param x,y,lab Column names for logFC, p-value, gene label (auto-detected
#'   if NULL).
#' @param custom_gmt Reserved for v0.2 (GMT path input).
#' @param p_method Y-axis transform: `"pi_eq2"` (default), `"pi_eq1"`,
#'   `"raw_p"`, `"adj_p"`.
#' @param p_adjust Method: `"BH"` (default), `"bonferroni"`, `"qvalue"`, `"IHW"`.
#' @param p_threshold,logfc_threshold Volcano significance cutoffs.
#' @param enrich_mode `c("fgsea", "ora")` - either or both.
#' @param enrich_padj Pathway significance cutoff (default 0.05).
#' @param dedup List with `method` (`"jaccard"`, `"collapse_fgsea"`, or
#'   `"both"`), `cutoff`, and `scope`. See [ev_collapse()].
#' @param ring List with `max_terms`, `order_by`, `magnitude`, `color`.
#' @param volcano List with `label_n`, `label_by`.
#' @param facet List with `nrow`, `ncol` for patchwork outer layout.
#' @param theme Output of `ev_theme()`.
#' @return Object of class `c("enrichVolcano", "patchwork")` with
#'   `attr(., "ev_data")` and `attr(., "ev_call")`.
#' @export
enrich_volcano <- function(data, contrast,
                           species = "human",
                           databases = c("hallmark", "reactome", "go_bp"),
                           x = NULL, y = NULL, lab = NULL,
                           custom_gmt = NULL,
                           p_method = c("pi_eq2", "pi_eq1", "raw_p", "adj_p"),
                           p_adjust = "BH",
                           p_threshold = 0.05,
                           logfc_threshold = log2(1.5),
                           enrich_mode = c("fgsea", "ora"),
                           enrich_padj = 0.05,
                           dedup = list(method = "jaccard", cutoff = 0.5,
                                        scope = "within_db"),
                           ring = list(max_terms = 10, order_by = "padj",
                                       magnitude = "neg_log_padj",
                                       color = "nes"),
                           volcano = list(label_n = 10, label_by = NULL),
                           facet = list(nrow = NULL, ncol = NULL),
                           theme = ev_theme()) {
  call <- match.call()
  p_method <- match.arg(p_method)

  validated <- ev_validate(data, x = x, y = y, lab = lab)

  with_pi <- if (p_method %in% c("pi_eq2", "pi_eq1")) {
    pi_score(validated, variant = sub("pi_", "", p_method))
  } else {
    validated
  }

  with_padj <- adjust_p(with_pi, method = p_adjust)

  # fgsea is ranked by a signed statistic (sign(logFC) * -log10(P)), NOT the
  # pi-score, so up/down genes sit at opposite ends of the ranked list.
  rank_by <- "signed_p"

  enrich <- ev_enrich(
    with_padj, contrast = contrast, databases = databases,
    species = species, enrich_mode = enrich_mode,
    rank_by = rank_by
  )

  if (nrow(enrich) > 0) {
    enrich <- enrich[enrich$padj < enrich_padj, , drop = FALSE]
  }

  dedup_res <- if (nrow(enrich) > 0) {
    ev_collapse(
      enrich,
      method = dedup$method, cutoff = dedup$cutoff,
      scope = dedup$scope
    )
  } else {
    enrich$dedup_kept <- logical(0)
    enrich
  }

  composites <- stats::setNames(
    lapply(contrast, function(ctr) {
      volc_sub <- with_padj[with_padj$contrast == ctr, , drop = FALSE]
      enr_sub <- dedup_res[dedup_res$contrast == ctr, , drop = FALSE]
      if ("dedup_kept" %in% colnames(enr_sub)) {
        enr_sub <- enr_sub[enr_sub$dedup_kept, , drop = FALSE]
      }
      if (nrow(enr_sub) > ring$max_terms) {
        enr_sub <- enr_sub[order(enr_sub$padj), , drop = FALSE]
        enr_sub <- utils::head(enr_sub, ring$max_terms)
      }
      ev_volcano_ring(volc_sub, enr_sub, title = ctr,
                      p_threshold = p_threshold, theme = theme)
    }),
    contrast
  )

  data_attr <- list(
    validated_input = validated,
    pi_scores = with_pi,
    enrichment = enrich,
    dedup_result = dedup_res
  )

  p <- if (length(composites) == 1) {
    composites[[1]]
  } else {
    patchwork::wrap_plots(composites, nrow = facet$nrow, ncol = facet$ncol)
  }
  if (!inherits(p, "enrichVolcano")) {
    class(p) <- c("enrichVolcano", class(p))
  }
  attr(p, "ev_data") <- data_attr
  attr(p, "ev_call") <- call
  p
}
