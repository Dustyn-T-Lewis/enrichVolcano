#' List the term labels a figure will draw, to review or edit
#'
#' Picks terms the way [plot_volcano_ring()] does, for each contrast, and lists
#' each term once with both shipped styles from [clean_label()]. Edit the
#' `label` column, in R or in the CSV that `write_labels()` writes, and pass
#' the table back through `labels` of [plot_volcano_ring()] or
#' [plot_scatter()]. A line break is written as the two characters `\n`.
#'
#' @inheritParams plot_volcano_ring
#' @param contrast Contrasts to list; `NULL` (default) lists every contrast.
#' @param n_terms Most unique terms per contrast, as on the ring. `Inf` lists
#'   every significant term.
#' @return A data frame with one row per term: `database`, `term`,
#'   `contrasts` (the contrasts that draw it, joined by `;`), `clean`, `short`
#'   and `label`, a copy of `short` to edit.
#' @export
#' @examples
#' ex <- as_enrichment(read.csv(system.file("extdata", "examples", "yvo_fgsea.csv.gz",
#'   package = "enrichVolcano"
#' )))
#' head(list_labels(ex))
list_labels <- function(enrichment, contrast = NULL, databases = NULL, collapse = TRUE,
                        term_threshold = 0.05, n_terms = 12) {
  check_enrichment(enrichment)
  res <- enrichment@results
  contrasts <- contrast %||% unique(res$contrast)
  picked <- do.call(rbind, lapply(contrasts, function(ct) {
    rows <- filter_view(contrast_rows(res, ct), databases, collapse)
    ring_terms(rows, term_threshold, n_terms, NULL)[, c("contrast", "database", "term", "padj")]
  }))
  picked <- picked[order(picked$padj), , drop = FALSE]
  terms <- unique(picked$term)
  flat <- function(style) gsub("\n", "\\n", clean_label(terms, width = 10000, style = style), fixed = TRUE)
  short <- flat("short")
  data.frame(
    database = picked$database[match(terms, picked$term)],
    term = terms,
    contrasts = vapply(terms, function(t) {
      paste(intersect(contrasts, picked$contrast[picked$term == t]), collapse = ";")
    }, character(1), USE.NAMES = FALSE),
    clean = flat("clean"),
    short = short,
    label = short
  )
}

#' @rdname list_labels
#' @param file Path of the CSV to write.
#' @param ... Arguments passed on to `list_labels()`.
#' @return `write_labels()`: `file`, invisibly.
#' @export
write_labels <- function(enrichment, file, ...) {
  utils::write.csv(list_labels(enrichment, ...), file, row.names = FALSE)
  invisible(file)
}
