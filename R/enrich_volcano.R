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
#' @param dedup List with `method`, `cutoff`, `scope`, `collapse_fgsea`.
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
                                        scope = "within_db",
                                        collapse_fgsea = TRUE),
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

  rank_by <- if (p_method == "pi_eq1") "pi_eq1" else "pi_eq2"

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

  volcano_plots <- stats::setNames(
    lapply(contrast, function(ctr) {
      ev_volcano(with_padj, contrast = ctr,
                 p_method = p_method,
                 p_threshold = p_threshold,
                 logfc_threshold = logfc_threshold,
                 label_n = volcano$label_n,
                 label_by = volcano$label_by,
                 theme = theme)
    }),
    contrast
  )

  ring_plots <- stats::setNames(
    lapply(contrast, function(ctr) {
      ring_plot(dedup_res, contrast = ctr,
              max_terms = ring$max_terms,
              order_by = ring$order_by,
              magnitude = ring$magnitude,
              color = ring$color,
              theme = theme)
    }),
    contrast
  )

  data_attr <- list(
    validated_input = validated,
    pi_scores = with_pi,
    enrichment = enrich,
    dedup_result = dedup_res
  )

  p <- ev_compose(volcano_plots, ring_plots,
                  nrow = facet$nrow, ncol = facet$ncol,
                  data = data_attr)
  attr(p, "ev_call") <- call
  p
}
