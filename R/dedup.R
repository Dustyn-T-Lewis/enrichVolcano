#' Flag redundant terms for display
#'
#' Related gene sets often reach significance together: in GO:BP or
#' Reactome, one biological signal can light up a dozen overlapping terms.
#' `dedup_terms()` picks one representative per cluster so a figure shows the
#' signal once. It changes what is drawn, not what was tested: no row is
#' dropped and no p-value is touched, because the multiple-testing correction
#' already covered every term.
#'
#' @section How terms are compared:
#' Within each contrast and database, significant terms are walked from the
#' smallest `padj` up. A term is `"redundant"` when its gene set is at least
#' `cutoff` similar to a representative already kept, and is recorded as
#' merged into the most similar one; otherwise it becomes a representative.
#' Similarity uses the full gene sets, not leading edges.
#'
#' * `"combined"` (default, cutoff 0.375) is the EnrichmentMap coefficient:
#'   the mean of the Jaccard index and the overlap coefficient
#'   (Merico et al. 2010; Reimand et al. 2019). The overlap coefficient
#'   catches a small set nested inside a large one, which Jaccard misses.
#' * `"jaccard"` (cutoff 0.5) uses shared genes over all genes in either set.
#'
#' @section `method = "collapse_pathways"`:
#' Runs [fgsea::collapsePathways()] within each contrast and database: a term
#' is redundant when it is no longer enriched once the genes of a more
#' significant term are conditioned on. This is a statistical criterion, not
#' an overlap rule, so it needs the ranking each contrast was tested on and
#' applies only to ranked GSEA results (fgsea, clusterProfiler). It permutes,
#' so call [set.seed()] first for reproducible flags. The p-values it computes
#' decide redundancy only; the reported `padj` stay those of the original run.
#'
#' @param enrichment An enrichment object from [as_enrichment()] or [run_enrichment()].
#' @param gene_sets A named list of character vectors, one per term, holding
#'   the full gene sets that were tested. Terms without a set are kept.
#' @param method `"enrichmentmap"` (gene-set overlap) or `"collapse_pathways"`
#'   (conditional enrichment).
#' @param similarity For `"enrichmentmap"`: `"combined"` or `"jaccard"`.
#' @param cutoff For `"enrichmentmap"`: similarity at or above which a term is
#'   redundant. `NULL` takes 0.375 for `"combined"` and 0.5 for `"jaccard"`.
#' @param term_threshold Only terms with `padj` below this are compared; the rest
#'   are left unflagged (`NA`). For `"collapse_pathways"` it is also the
#'   conditional p-value threshold.
#' @param stats For `"collapse_pathways"`: a named list with one ranking (named
#'   numeric vector of gene statistics) per contrast.
#'
#' @return `enrichment`, with `dedup_status` (`"kept"`, `"redundant"` or `NA`),
#'   `merged_into` and `similarity` columns in its results, and the settings
#'   stored in `metadata$dedup`.
#' @references
#' Merico D, Isserlin R, Stueker O, Emili A, Bader GD (2010). Enrichment Map: a
#' network-based method for gene-set enrichment visualization and
#' interpretation. PLoS ONE 5(11):e13984.
#'
#' Reimand J, Isserlin R, Voisin V, et al. (2019). Pathway enrichment analysis
#' and visualization of omics data using g:Profiler, GSEA, Cytoscape and
#' EnrichmentMap. Nature Protocols 14:482-517.
#' @export
dedup_terms <- function(enrichment, gene_sets, method = c("enrichmentmap", "collapse_pathways"),
                        similarity = c("combined", "jaccard"), cutoff = NULL,
                        term_threshold = 0.05, stats = NULL) {
  check_enrichment(enrichment)
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
    ev_abort("{.arg gene_sets} must be a named list of gene vectors, one per term.",
      class = "enrichVolcano_input_error"
    )
  }
  method <- rlang::arg_match(method)
  res <- enrichment@results
  res[intersect(c("dedup_status", "merged_into", "similarity", "overlap_jaccard"), names(res))] <- NULL

  if (method == "enrichmentmap") {
    similarity <- rlang::arg_match(similarity)
    cutoff <- cutoff %||% c(combined = 0.375, jaccard = 0.5)[[similarity]]
    if (!is.numeric(cutoff) || length(cutoff) != 1 || is.na(cutoff) || cutoff <= 0 || cutoff > 1) {
      ev_abort("{.arg cutoff} must be a single number in (0, 1].", class = "enrichVolcano_param_error")
    }
    sim <- switch(similarity,
      combined = combined_similarity,
      jaccard = jaccard
    )
    flags <- flag_redundant(res, gene_sets, sim, cutoff, term_threshold)
  } else {
    check_collapse_inputs(enrichment, stats)
    similarity <- NULL
    cutoff <- NULL
    flags <- flag_collapsed(res, gene_sets, stats, term_threshold)
  }

  enrichment@results <- cbind(res, flags)
  enrichment@metadata$dedup <- list(
    method = method, similarity = similarity, cutoff = cutoff, term_threshold = term_threshold
  )
  enrichment
}

check_collapse_inputs <- function(x, stats) {
  if (!x@metadata$enrichment_test %in% c("fgsea", "gseaResult")) {
    ev_abort(
      c(
        "{.val collapse_pathways} needs ranked GSEA results, not {.val {x@metadata$enrichment_test}}.",
        i = "Use {.code method = \"enrichmentmap\"} for these."
      ),
      class = "enrichVolcano_input_error"
    )
  }
  if (is.null(stats)) {
    ev_abort(
      "{.val collapse_pathways} needs {.arg stats}: the ranking each contrast was tested on.",
      class = "enrichVolcano_param_error"
    )
  }
  missing <- setdiff(unique(x@results$contrast), names(stats))
  if (length(missing) > 0) {
    ev_abort("{.arg stats} has no ranking for contrast{?s} {.val {missing}}.",
      class = "enrichVolcano_input_error"
    )
  }
  rlang::check_installed(c("fgsea", "data.table"), reason = "for `method = \"collapse_pathways\"`.")
}

flag_collapsed <- function(res, gene_sets, stats, term_threshold) {
  status <- rep(NA_character_, nrow(res))
  merged_into <- rep(NA_character_, nrow(res))

  sig <- significant_rows(res, gene_sets, term_threshold)
  status[sig] <- "kept"
  sig <- sig[res$term[sig] %in% names(gene_sets)]
  for (rows in split(sig, paste(res$contrast[sig], res$database[sig], sep = "\r"))) {
    rows <- rows[order(res$p[rows])]
    fg <- data.table::data.table(
      pathway = res$term[rows], ES = res$score[rows],
      pval = res$p[rows], padj = res$padj[rows]
    )
    collapsed <- fgsea::collapsePathways(
      fg, gene_sets, stats[[res$contrast[rows[1]]]],
      pval.threshold = term_threshold
    )
    parent <- collapsed$parentPathways[res$term[rows]]
    status[rows] <- ifelse(is.na(parent), "kept", "redundant")
    merged_into[rows] <- unname(parent)
  }
  data.frame(dedup_status = status, merged_into = merged_into, similarity = NA_real_)
}

flag_redundant <- function(res, gene_sets, sim, cutoff, term_threshold) {
  status <- rep(NA_character_, nrow(res))
  merged_into <- rep(NA_character_, nrow(res))
  similarity <- rep(NA_real_, nrow(res))

  sig <- significant_rows(res, gene_sets, term_threshold)

  for (rows in split(sig, paste(res$contrast[sig], res$database[sig], sep = "\r"))) {
    kept <- integer(0)
    for (i in rows[order(res$padj[rows])]) {
      genes <- gene_sets[[res$term[i]]]
      s <- vapply(kept, function(k) sim(genes, gene_sets[[res$term[k]]]), numeric(1))
      if (!is.null(genes) && length(s) > 0 && max(s) >= cutoff) {
        best <- which.max(s)
        status[i] <- "redundant"
        merged_into[i] <- res$term[kept[best]]
        similarity[i] <- s[best]
      } else {
        status[i] <- "kept"
        if (!is.null(genes)) kept <- c(kept, i)
      }
    }
  }
  data.frame(dedup_status = status, merged_into = merged_into, similarity = similarity)
}

significant_rows <- function(res, gene_sets, term_threshold) {
  sig <- which(!is.na(res$padj) & res$padj < term_threshold)
  no_set <- setdiff(res$term[sig], names(gene_sets))
  if (length(no_set) > 0) {
    ev_inform("{length(no_set)} term{?s} had no gene set and {?was/were} kept.",
      class = "enrichVolcano_dedup_missing_sets"
    )
  }
  sig
}

jaccard <- function(a, b) {
  shared <- length(intersect(a, b))
  if (shared == 0) {
    return(0)
  }
  shared / length(union(a, b))
}

overlap_coefficient <- function(a, b) {
  shared <- length(intersect(a, b))
  if (shared == 0) {
    return(0)
  }
  shared / min(length(unique(a)), length(unique(b)))
}

combined_similarity <- function(a, b) {
  0.5 * jaccard(a, b) + 0.5 * overlap_coefficient(a, b)
}
