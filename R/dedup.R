#' Flag redundant terms for display
#'
#' Related gene sets often reach significance together: in GO:BP or
#' Reactome, one biological signal can light up a dozen overlapping terms.
#' `dedup()` picks one representative per cluster so a figure shows the
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
#' @param x An [enrichment] object.
#' @param gene_sets A named list of character vectors, one per term, holding
#'   the full gene sets that were tested. Terms without a set are kept.
#' @param method `"enrichmentmap"`.
#' @param similarity `"combined"` or `"jaccard"`.
#' @param cutoff Similarity at or above which a term is redundant. `NULL`
#'   takes 0.375 for `"combined"` and 0.5 for `"jaccard"`.
#' @param p_threshold Only terms with `padj` below this are compared; the rest
#'   are left unflagged (`NA`).
#'
#' @return `x`, with `dedup_status` (`"kept"`, `"redundant"` or `NA`),
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
dedup <- function(x, gene_sets, method = "enrichmentmap",
                  similarity = c("combined", "jaccard"), cutoff = NULL,
                  p_threshold = 0.05) {
  if (!S7::S7_inherits(x, enrichment)) {
    ev_abort("{.arg x} must be an {.cls enrichment}; build one with {.fn as_enrichment}.",
      class = "enrichVolcano_input_error"
    )
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
    ev_abort("{.arg gene_sets} must be a named list of gene vectors, one per term.",
      class = "enrichVolcano_input_error"
    )
  }
  method <- rlang::arg_match(method, "enrichmentmap")
  similarity <- rlang::arg_match(similarity)
  cutoff <- cutoff %||% c(combined = 0.375, jaccard = 0.5)[[similarity]]
  if (!is.numeric(cutoff) || length(cutoff) != 1 || is.na(cutoff) || cutoff <= 0 || cutoff > 1) {
    ev_abort("{.arg cutoff} must be a single number in (0, 1].", class = "enrichVolcano_param_error")
  }

  res <- x@results
  res[intersect(c("dedup_status", "merged_into", "similarity", "overlap_jaccard"), names(res))] <- NULL
  sim <- switch(similarity,
    combined = combined_similarity,
    jaccard = jaccard
  )
  x@results <- cbind(res, flag_redundant(res, gene_sets, sim, cutoff, p_threshold))
  x@metadata$dedup <- list(
    method = method, similarity = similarity, cutoff = cutoff, p_threshold = p_threshold
  )
  x
}

flag_redundant <- function(res, gene_sets, sim, cutoff, p_threshold) {
  status <- rep(NA_character_, nrow(res))
  merged_into <- rep(NA_character_, nrow(res))
  similarity <- rep(NA_real_, nrow(res))

  sig <- which(!is.na(res$padj) & res$padj < p_threshold)
  no_set <- setdiff(res$term[sig], names(gene_sets))
  if (length(no_set) > 0) {
    ev_inform("{length(no_set)} term{?s} had no gene set and {?was/were} kept.",
      class = "enrichVolcano_dedup_missing_sets"
    )
  }

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
