#' Run gene-set enrichment on a study
#'
#' Tests each contrast of a study against gene-set collections and returns one
#' [enrichment] object per test, ready for the plots.
#'
#' @section fgsea:
#' Proteins are ranked by the study's `rank` column (see [as_da()]). Where
#' several proteins share a gene symbol, the most abundant one represents the
#' gene (the table's average abundance, else the matrix row mean); without any
#' abundance, the first accession alphabetically. The choice never looks at
#' the results. Each collection is tested and corrected separately. fgsea
#' permutes, so call [set.seed()] first for reproducible p-values.
#'
#' @param study Output of [read_study()], or of [as_da()] for fgsea alone.
#' @param gene_sets Output of [load_gene_sets()], a named list of collections
#'   of gene-symbol vectors, or one named list of sets.
#' @param tests `"fgsea"`.
#' @param min_size,max_size Sets need this many genes present in the data.
#' @return A named list with one [enrichment] object per test. Its metadata
#'   records the ranking statistic and the gene-set versions.
#' @export
run_enrichment <- function(study, gene_sets, tests = "fgsea", min_size = 15, max_size = 500) {
  study <- as_study(study)
  tests <- rlang::arg_match(tests, "fgsea", multiple = TRUE)
  collections <- as_collections(gene_sets)
  out <- list()
  if ("fgsea" %in% tests) out$fgsea <- run_fgsea(study, collections, min_size, max_size)
  out
}

as_study <- function(x) {
  if (inherits(x, "enrichVolcano_study")) {
    return(x)
  }
  if (inherits(x, "enrichVolcano_da")) {
    return(structure(list(da = x), class = c("enrichVolcano_study", "list")))
  }
  ev_abort(
    "{.arg study} must come from {.fn read_study} or {.fn as_da}, not {.cls {class(x)[1]}}.",
    class = "enrichVolcano_input_error"
  )
}

as_collections <- function(gene_sets) {
  is_set_list <- function(x) {
    is.list(x) && length(x) > 0 && !is.null(names(x)) && all(vapply(x, is.character, logical(1)))
  }
  if (is_set_list(gene_sets)) {
    return(structure(list(`gene sets` = gene_sets), versions = attr(gene_sets, "versions")))
  }
  if (is.list(gene_sets) && !is.null(names(gene_sets)) && all(vapply(gene_sets, is_set_list, logical(1)))) {
    return(gene_sets)
  }
  ev_abort(
    "{.arg gene_sets} must be a named list of gene-symbol vectors, or a named list of such collections.",
    class = "enrichVolcano_input_error"
  )
}

run_fgsea <- function(study, collections, min_size, max_size) {
  rlang::check_installed("fgsea", reason = "to run fgsea.")
  da <- fill_abundance(study$da, study$matrix)
  da <- da[!is.na(da$gene) & !is.na(da$rank), , drop = FALSE]
  per_contrast <- lapply(split(da, factor(da$contrast, unique(da$contrast))), function(d) {
    d <- one_per_gene(d)
    ranks <- stats::setNames(d$rank, d$gene)
    do.call(rbind, lapply(names(collections), function(db) {
      res <- as.data.frame(fgsea::fgsea(collections[[db]], ranks, minSize = min_size, maxSize = max_size))
      res$database <- rep(db, nrow(res))
      res
    }))
  })
  if (sum(vapply(per_contrast, nrow, integer(1))) == 0) {
    ev_abort(
      "There is no gene set with {min_size} to {max_size} genes present in the data.",
      class = "enrichVolcano_input_error"
    )
  }
  x <- suppressMessages(as_enrichment(per_contrast, enrichment_test = "fgsea"))
  x@metadata$ranking <- unique(da$rank_stat)
  x@metadata$gene_sets <- attr(collections, "versions")
  x
}

one_per_gene <- function(d) {
  d <- d[order(d$gene, -d$abundance, d$protein, na.last = TRUE), , drop = FALSE]
  d[!duplicated(d$gene), , drop = FALSE]
}

fill_abundance <- function(da, matrix) {
  if (is.null(matrix)) {
    return(da)
  }
  means <- rowMeans(matrix)
  missing <- is.na(da$abundance) & da$protein %in% names(means)
  da$abundance[missing] <- means[da$protein[missing]]
  da
}
