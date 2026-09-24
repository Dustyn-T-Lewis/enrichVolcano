#' Run gene-set enrichment on a study
#'
#' Tests each contrast of a study against gene-set collections and returns one
#' enrichment object per test, the same object [as_enrichment()] builds.
#'
#' @section fgsea:
#' Proteins are ranked by the study's `rank` column (see [as_da()]). Where
#' several proteins share a gene symbol, the most abundant one represents the
#' gene (the table's average abundance, else the matrix row mean); without any
#' abundance, the first accession alphabetically. The choice never looks at
#' the results. Each collection is tested and corrected separately. fgsea
#' permutes, so call [set.seed()] first for reproducible p-values.
#'
#' @section camera and fry:
#' Both refit limma on the study's matrix, which must have no missing values
#' (use an imputed or limpa-quantified matrix), reduced to one protein per gene
#' by the same rule. The design is `~ 0 + group`, plus any `covariates`, plus
#' `subject` when `subject_effect = "fixed"`; contrasts come from the
#' `contrasts` sheet. With `subject_effect = "block"` and a `subject` column,
#' samples from one subject are treated as correlated
#' (`limma::duplicateCorrelation()`). fry uses that blocking directly. camera
#' cannot block, so for a blocked design the package runs
#' `limma::cameraPR()` on the moderated t of the blocked fit instead, and the
#' result records `cameraPR` as its test. Precision weights are used when the
#' study has them. camera assumes genes within a set correlate at
#' `inter_gene_cor` (limma's default, 0.01); `NA` estimates it from the data,
#' which can change results substantially and needs an unblocked design.
#'
#' @param study Output of [read_study()], or of [as_da()] for fgsea alone.
#' @param gene_sets Output of [load_gene_sets()], a named list of collections
#'   of gene-symbol vectors, or one named list of sets.
#' @param tests Any of `"fgsea"`, `"camera"`, `"fry"`. A study without a
#'   matrix runs fgsea only.
#' @param subject_effect `"block"` or `"fixed"`: how a `subject` column enters
#'   the camera and fry model.
#' @param covariates Sample columns added to the camera and fry design.
#' @param inter_gene_cor camera's inter-gene correlation: a number, or `NA` to
#'   estimate it (unblocked designs only).
#' @param min_size,max_size Sets need this many genes present in the data.
#' @return A named list with one enrichment object ([as_enrichment()]) per test. Its metadata
#'   records the ranking statistic and the gene-set versions.
#' @export
run_enrichment <- function(study, gene_sets, tests = c("fgsea", "camera", "fry"),
                           subject_effect = c("block", "fixed"), covariates = NULL,
                           inter_gene_cor = 0.01, min_size = 15, max_size = 500) {
  study <- as_study(study)
  tests_given <- !missing(tests)
  tests <- rlang::arg_match(tests, c("fgsea", "camera", "fry"), multiple = TRUE)
  subject_effect <- rlang::arg_match(subject_effect)
  collections <- as_collections(gene_sets)
  limma_tests <- intersect(c("camera", "fry"), tests)
  if (is.null(study$matrix) && length(limma_tests) > 0) {
    if (tests_given) {
      ev_abort(
        "camera and fry need a matrix, samples and contrasts; see {.fn read_study}.",
        class = "enrichVolcano_input_error"
      )
    }
    ev_inform("This study has no matrix, so only fgsea runs.", class = "enrichVolcano_fgsea_only")
    limma_tests <- character(0)
  }
  out <- list()
  if ("fgsea" %in% tests) out$fgsea <- run_fgsea(study, collections, min_size, max_size)
  if (length(limma_tests) > 0) {
    model <- limma_model(study, subject_effect, covariates)
    model$inter_gene_cor <- inter_gene_cor
    if ("camera" %in% limma_tests && is.na(inter_gene_cor) && !is.null(model$block)) {
      ev_abort(
        c(
          "A blocked design runs camera as {.fn limma::cameraPR}, which cannot estimate {.arg inter_gene_cor}.",
          i = "Give a number, or use {.code subject_effect = \"fixed\"}."
        ),
        class = "enrichVolcano_param_error"
      )
    }
    for (test in limma_tests) out[[test]] <- run_limma_test(test, model, collections, min_size, max_size)
  }
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
  ties <- 0
  per_contrast <- lapply(split(da, factor(da$contrast, unique(da$contrast))), function(d) {
    d <- one_per_gene(d)
    ranks <- stats::setNames(d$rank, d$gene)
    tied <- duplicated(ranks) | duplicated(ranks, fromLast = TRUE)
    if (sum(tied) > ties) ties <<- sum(tied)
    do.call(rbind, lapply(names(collections), function(db) {
      res <- withCallingHandlers(
        as.data.frame(fgsea::fgsea(collections[[db]], ranks, minSize = min_size, maxSize = max_size)),
        warning = function(w) {
          if (grepl("ties in the preranked stats", conditionMessage(w))) invokeRestart("muffleWarning")
        }
      )
      res$database <- rep(db, nrow(res))
      res
    }))
  })
  if (ties > 0) {
    ev_inform(
      "Up to {ties} of {length(unique(da$gene))} ranked genes tie; fgsea orders tied genes arbitrarily.",
      class = "enrichVolcano_rank_ties"
    )
  }
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
  means <- rowMeans(matrix, na.rm = TRUE)
  missing <- is.na(da$abundance) & da$protein %in% names(means)
  da$abundance[missing] <- means[da$protein[missing]]
  da
}

limma_model <- function(study, subject_effect, covariates) {
  rlang::check_installed("limma", reason = "to run camera and fry.")
  n_missing <- sum(is.na(study$matrix))
  if (n_missing > 0) {
    ev_abort(
      c(
        "The matrix has {n_missing} missing value{?s}.",
        i = "camera and fry need a complete matrix: use an imputed or limpa-quantified one."
      ),
      class = "enrichVolcano_data_error"
    )
  }
  genes <- fill_abundance(study$da, study$matrix)
  genes <- genes[!duplicated(genes$protein) & !is.na(genes$gene) & genes$protein %in% rownames(study$matrix), ]
  genes <- one_per_gene(genes)
  m <- study$matrix[genes$protein, , drop = FALSE]
  rownames(m) <- genes$gene
  w <- study$weights
  if (!is.null(w)) {
    w <- w[genes$protein, , drop = FALSE]
    rownames(w) <- genes$gene
  }
  samples <- study$samples
  design <- design_matrix(samples, subject_effect, covariates)
  cm <- limma::makeContrasts(contrasts = study$contrasts$expression, levels = design)
  colnames(cm) <- study$contrasts$name

  model <- list(m = m, w = w, design = design, cm = cm, block = NULL, correlation = NULL)
  if (subject_effect == "block") {
    if (is.null(samples$subject)) {
      ev_inform("No {.field subject} column, so samples are treated as independent.",
        class = "enrichVolcano_no_blocking"
      )
    } else {
      model$block <- samples$subject
      dupcor <- limma::duplicateCorrelation(m, design, block = model$block, weights = w)
      model$correlation <- dupcor$consensus.correlation
      fit <- limma::lmFit(m, design, block = model$block, correlation = model$correlation, weights = w)
      model$t <- limma::eBayes(limma::contrasts.fit(fit, cm))$t
    }
  }
  model
}

design_matrix <- function(samples, subject_effect, covariates) {
  require_columns(samples, covariates)
  if (subject_effect == "fixed") require_columns(samples, "subject")
  group <- factor(samples$group)
  data <- data.frame(group = group, samples[covariates], check.names = FALSE)
  terms <- c("0", "group", covariates)
  if (subject_effect == "fixed") {
    data$subject <- factor(samples$subject)
    terms <- c(terms, "subject")
  }
  design <- stats::model.matrix(stats::reformulate(terms), data)
  colnames(design)[seq_len(nlevels(group))] <- levels(group)
  design
}

run_limma_test <- function(test, model, collections, min_size, max_size) {
  blocked_camera <- test == "camera" && !is.null(model$block)
  per_contrast <- lapply(colnames(model$cm), function(contrast) {
    do.call(rbind, lapply(names(collections), function(db) {
      index <- limma::ids2indices(collections[[db]], rownames(model$m))
      index <- index[lengths(index) >= min_size & lengths(index) <= max_size]
      if (length(index) == 0) {
        return(NULL)
      }
      res <- if (test == "fry") {
        limma::fry(model$m, index, model$design, model$cm[, contrast],
          block = model$block, correlation = model$correlation, weights = model$w, sort = "none"
        )
      } else if (blocked_camera) {
        limma::cameraPR(model$t[, contrast], index, inter.gene.cor = model$inter_gene_cor, sort = FALSE)
      } else {
        limma::camera(model$m, index, model$design, model$cm[, contrast],
          weights = model$w, inter.gene.cor = model$inter_gene_cor, sort = FALSE
        )
      }
      res$term <- rownames(res)
      res$database <- db
      rownames(res) <- NULL
      res
    }))
  })
  names(per_contrast) <- colnames(model$cm)
  if (all(vapply(per_contrast, is.null, logical(1)))) {
    ev_abort(
      "There is no gene set with {min_size} to {max_size} genes present in the data.",
      class = "enrichVolcano_input_error"
    )
  }
  name <- if (blocked_camera) "cameraPR" else test
  x <- as_enrichment(per_contrast, enrichment_test = name)
  x@metadata$gene_sets <- attr(collections, "versions")
  if (test == "camera") x@metadata$inter_gene_cor <- model$inter_gene_cor
  x@metadata$design <- list(
    columns = colnames(model$design), blocked = !is.null(model$block),
    correlation = model$correlation, weighted = !is.null(model$w)
  )
  x
}
