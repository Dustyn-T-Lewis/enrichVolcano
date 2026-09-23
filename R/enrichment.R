enrichment_tests <- c(
  "fgsea", "fry", "mroast", "camera", "cameraPR", "gseaResult", "ora", "custom"
)

results_columns <- c(
  contrast = "character", database = "character", term = "character",
  score = "numeric", p = "numeric", padj = "numeric", size = "numeric",
  direction = "character", leading_edge = "list"
)

#' Enrichment results in one validated object
#'
#' `enrichment` is the single input every plot in the package reads. Build it
#' with [as_enrichment()] rather than by hand; construct it directly only when
#' you already hold a table in exactly this shape.
#'
#' @param results A data frame with one row per term per contrast:
#'   * `contrast`, `database`, `term`: character. `database` may be `NA`.
#'   * `score`: numeric. NES for fgsea and clusterProfiler; signed
#'     \eqn{-\log_{10}}(FDR) for limma rotation/competitive tests and ORA.
#'   * `p`, `padj`: numeric in \[0, 1\]; `NA` allowed.
#'   * `size`: numeric set size.
#'   * `direction`: `"up"` or `"down"`, agreeing with the sign of `score`.
#'   * `leading_edge`: list of character vectors, empty when the test has none.
#'
#'   Extra columns (for example `dedup_status` or your own groupings) are kept
#'   and can be mapped by the plots.
#' @param metadata A list with `enrichment_test` (the method that produced the
#'   results), `score_type` (the axis and legend label for `score`), and
#'   `dedup` (`NULL`, or the settings `dedup()` used).
#'
#' @return An S7 object with properties `results` and `metadata`. Every
#'   construction and every edit is validated.
#' @export
enrichment <- S7::new_class(
  "enrichment",
  package = "enrichVolcano",
  properties = list(
    results = S7::new_property(S7::class_data.frame, default = quote(data.frame())),
    metadata = S7::new_property(S7::class_list, default = quote(list()))
  ),
  validator = function(self) {
    enrichment_problems(self@results, self@metadata)
  }
)

enrichment_problems <- function(results, metadata) {
  missing <- setdiff(names(results_columns), names(results))
  if (length(missing) > 0) {
    return(paste0("`results` is missing ", paste0("`", missing, "`", collapse = ", "), "."))
  }
  wrong_type <- Filter(
    function(col) !is_column_type(results[[col]], results_columns[[col]]),
    names(results_columns)
  )
  if (length(wrong_type) > 0) {
    return(vapply(wrong_type, function(col) {
      sprintf(
        "`%s` must be %s, not %s.",
        col, results_columns[[col]], class(results[[col]])[1]
      )
    }, character(1), USE.NAMES = FALSE))
  }

  problems <- c(
    if (any(!in_unit_interval(results$padj))) "`padj` must lie in [0, 1].",
    if (any(!in_unit_interval(results$p))) "`p` must lie in [0, 1].",
    direction_problem(results$direction, results$score),
    duplicate_problem(results),
    metadata_problems(metadata)
  )
  if (length(problems) > 0) problems
}

is_column_type <- function(x, type) {
  switch(type,
    character = is.character(x),
    numeric = is.numeric(x),
    list = is.list(x)
  )
}

in_unit_interval <- function(x) is.na(x) | (x >= 0 & x <= 1)

direction_problem <- function(direction, score) {
  if (!all(direction %in% c("up", "down"))) {
    return('`direction` must be "up" or "down".')
  }
  expected <- ifelse(score > 0, "up", "down")
  clash <- !is.na(score) & score != 0 & direction != expected
  if (any(clash)) {
    sprintf("`direction` disagrees with the sign of `score` in %d row(s).", sum(clash))
  }
}

duplicate_problem <- function(results) {
  key <- paste(results$contrast, results$database, results$term, sep = "\r")
  n_dup <- sum(duplicated(key))
  if (n_dup > 0) {
    sprintf("`results` has %d duplicate contrast / database / term row(s).", n_dup)
  }
}

metadata_problems <- function(metadata) {
  c(
    if (!isTRUE(metadata$enrichment_test %in% enrichment_tests)) {
      paste0(
        "`metadata$enrichment_test` must be one of ",
        paste0('"', enrichment_tests, '"', collapse = ", "), "."
      )
    },
    if (!rlang::is_string(metadata$score_type) || !nzchar(metadata$score_type)) {
      "`metadata$score_type` must be a non-empty string."
    }
  )
}

S7::method(print, enrichment) <- function(x, ...) {
  res <- x@results
  meta <- x@metadata
  contrasts <- unique(res$contrast)
  db <- table(ifelse(is.na(res$database), "unlabelled", res$database))
  dedup <- if (is.null(meta$dedup)) {
    "not deduplicated"
  } else {
    paste(meta$dedup$method, meta$dedup$similarity, meta$dedup$cutoff)
  }
  cli::cli_text("{.cls enrichment} {meta$enrichment_test}, score: {meta$score_type}")
  cli::cli_text("{length(contrasts)} contrast{?s}: {contrasts}")
  cli::cli_text("{nrow(res)} row{?s}: {paste0(names(db), ' (', db, ')')}")
  cli::cli_text("Dedup: {dedup}")
  invisible(x)
}
