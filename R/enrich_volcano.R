# Choose which enrichment rows the ring displays for one contrast.
ev_select_ring_terms <- function(enr, max_terms, show_databases = NULL,
                                  direction_balance = FALSE) {
  if (!is.null(show_databases)) {
    enr <- enr[enr$database %in% show_databases, , drop = FALSE]
  }
  if (nrow(enr) == 0) return(enr)
  enr <- enr[order(enr$padj), , drop = FALSE]
  if (!direction_balance) {
    return(utils::head(enr, max_terms))
  }
  up <- enr[!is.na(enr$NES) & enr$NES > 0, , drop = FALSE]
  dn <- enr[!is.na(enr$NES) & enr$NES < 0, , drop = FALSE]
  n_up <- min(nrow(up), ceiling(max_terms / 2))
  n_dn <- min(nrow(dn), max_terms - n_up)
  n_up <- min(nrow(up), max_terms - n_dn)  # let one side fill remainder
  out <- rbind(utils::head(up, n_up), utils::head(dn, n_dn))
  out[order(out$padj), , drop = FALSE]
}

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
#'   "fly", "yeast", "pig". UniProt accession mapping is available for `human`,
#'   `mouse`, and `rat`; other species are supported for enrichment but require
#'   gene-symbol input.
#' @param databases Character vector of registered DB names (see
#'   `list_databases()`), or a named list of pathway lists.
#' @param x,y,lab Column names for logFC, p-value, gene label (auto-detected
#'   if NULL).
#' @param uniprot Column name holding UniProt accessions (auto-detected by default).
#' @param p_method Y-axis transform: `"pi_eq2"` (default), `"pi_eq1"`,
#'   `"raw_p"`, `"adj_p"`.
#' @param rank_by Signed statistic fgsea ranks genes by. `"signed_p"`
#'   (default, `sign(logFC) * -log10(P)`), `"t"`, `"logFC"`, or any signed
#'   column. Unsigned pi-score columns are rejected.
#' @param p_adjust Method: `"BH"` (default), `"bonferroni"`, `"qvalue"`, `"IHW"`.
#' @param p_threshold,logfc_threshold Volcano significance cutoffs. A point is
#'   coloured up/down only when it clears both `p_threshold` and
#'   `abs(logFC) >= logfc_threshold`.
#' @param enrich_mode `c("fgsea", "ora")` - either or both.
#' @param enrich_padj Pathway significance cutoff (default 0.05).
#' @param dedup List with `method` (`"jaccard"`, `"collapse_fgsea"`, or
#'   `"both"`), `cutoff`, and `scope`. See [ev_collapse()].
#' @param ring List controlling the ring. Honoured fields: `max_terms` (cap on
#'   pathways drawn in the ring); `show_databases` (default NULL = all enriched
#'   databases; pass a character vector to restrict the ring to those named
#'   databases); `direction_balance` (default FALSE = global top-N by padj;
#'   TRUE = balanced selection of top up-NES + top down-NES terms; rows without
#'   a signed NES, such as ORA results, are excluded from balanced selection).
#'   Arc height
#'   is fixed to `-log10(padj)` and fill to NES by the canonical layout. For
#'   tunable `order_by`/`magnitude`/`color` encodings use the standalone
#'   [ring_plot()].
#' @param facet List with `nrow`, `ncol` for patchwork outer layout.
#' @param label_mode,label_n,label_rank_by,label_genes Volcano point labeling,
#'   passed to the shared selector; default `label_mode = "none"`.
#' @param disc_color Optional contrast-tint for the ring's central disc.
#' @param theme Output of `ev_theme()`.
#' @return Object of class `c("enrichVolcano", "patchwork")` with
#'   `attr(., "ev_data")` and `attr(., "ev_call")`.
#' @export
enrich_volcano <- function(data, contrast,
                           species = "human",
                           databases = c("hallmark", "reactome", "go_bp"),
                           x = NULL, y = NULL, lab = NULL,
                           uniprot = NULL,
                           p_method = c("pi_eq2", "pi_eq1", "raw_p", "adj_p"),
                           p_adjust = "BH",
                           p_threshold = 0.05,
                           logfc_threshold = log2(1.5),
                           enrich_mode = c("fgsea", "ora"),
                           enrich_padj = 0.05,
                           rank_by = "signed_p",
                           dedup = list(method = "jaccard", cutoff = 0.5,
                                        scope = "within_db"),
                           ring = list(max_terms = 10),
                           facet = list(nrow = NULL, ncol = NULL),
                           disc_color = NULL,
                           label_mode = c("none", "top_per_direction",
                                          "top_total", "all_significant",
                                          "explicit"),
                           label_n = 10,
                           label_rank_by = c("significance", "logfc"),
                           label_genes = NULL,
                           theme = ev_theme()) {
  call <- match.call()
  p_method <- match.arg(p_method)
  label_mode <- match.arg(label_mode)
  label_rank_by <- match.arg(label_rank_by)

  validated <- ev_validate(data, x = x, y = y, lab = lab,
                           uniprot = uniprot, species = species)

  with_pi <- if (p_method %in% c("pi_eq2", "pi_eq1")) {
    pi_score(validated, variant = sub("pi_", "", p_method))
  } else {
    validated
  }

  with_padj <- adjust_p(with_pi, method = p_adjust)

  # fgsea needs a SIGNED ranking statistic (default sign(logFC) * -log10(P)).
  # ev_enrich() rejects unsigned pi-score columns.
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
      enr_sub <- ev_select_ring_terms(
        enr_sub,
        max_terms         = ring$max_terms,
        show_databases    = ring$show_databases,
        direction_balance = isTRUE(ring$direction_balance)
      )
      ev_volcano_ring(volc_sub, enr_sub, title = ctr,
                      p_threshold = p_threshold,
                      logfc_threshold = logfc_threshold,
                      disc_color = disc_color, theme = theme,
                      label_mode = label_mode, label_n = label_n,
                      label_rank_by = label_rank_by, label_genes = label_genes,
                      p_method = p_method)
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
