#' Write enrichment results to a CSV file
#'
#' Writes one row per term per contrast, with an `enrichment_test` column and
#' leading edges joined by `;`. [as_enrichment()] reads the file back with
#' `enrichment_test = "custom"` and `leading_edge = "leading_edge"`.
#'
#' @param x An enrichment object, or the list of them that [run_enrichment()]
#'   returns, which is written as one long table.
#' @param file Path of the CSV to write.
#' @return `file`, invisibly.
#' @export
#' @examples
#' ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
#'   package = "enrichVolcano"
#' )))
#' write_table(ex, tempfile(fileext = ".csv"))
write_table <- function(x, file) {
  objects <- if (S7::S7_inherits(x, enrichment)) list(x) else x
  if (!is.list(objects) || length(objects) == 0 ||
    !all(vapply(objects, S7::S7_inherits, logical(1), enrichment))) {
    ev_abort(
      "{.arg x} must be an enrichment object or a list of them, as {.fn run_enrichment} returns.",
      class = "enrichVolcano_input_error"
    )
  }
  tables <- lapply(objects, function(e) {
    res <- e@results
    res$leading_edge <- vapply(res$leading_edge, paste, character(1), collapse = ";")
    cbind(enrichment_test = e@metadata$enrichment_test, res)
  })
  utils::write.csv(bind_rows_fill(tables), file, row.names = FALSE)
  invisible(file)
}
