#' Read per-contrast CSVs into a single validated long table
#'
#' Each CSV is one contrast; the file name (without extension) becomes the
#' contrast label. Required columns: a UniProt accession column, `logFC`, and a
#' p-value. An adjusted p-value column is optional; when absent it is computed
#' with [adjust_p()]. Accessions are mapped to gene symbols via the bundled
#' per-species table (see [ev_idmap_report()] for the mapping summary).
#'
#' @param path A directory (all matching files are read) or a character vector
#'   of file paths.
#' @param species One of `"human"`, `"mouse"`, `"rat"`.
#' @param p_adjust Adjustment method passed to [adjust_p()] when no adjusted
#'   column is present (default `"BH"`).
#' @param uniprot,x,y Optional column-name overrides for the accession, logFC,
#'   and p-value columns (passed to [ev_validate()]).
#' @param padj Optional name of an existing adjusted-p column to use as
#'   `adj.P.Val`.
#' @param pattern File-matching regex when `path` is a directory.
#' @return A validated long tibble with `uniprot`, `symbol`, `gene`, `contrast`,
#'   `logFC`, `P.Value`, `adj.P.Val`, and an `ev_idmap_report` attribute.
#' @export
ev_read_contrasts <- function(path, species, p_adjust = "BH",
                              uniprot = NULL, x = NULL, y = NULL,
                              padj = NULL, pattern = "\\.csv$") {
  files <- if (length(path) == 1 && dir.exists(path)) {
    list.files(path, pattern = pattern, full.names = TRUE)
  } else {
    path
  }
  if (length(files) == 0) {
    ev_abort("No contrast files found at {.path {path}}.",
             class = "ev_no_contrast_files")
  }
  labels <- tools::file_path_sans_ext(basename(files))
  if (anyDuplicated(labels)) {
    dup <- unique(labels[duplicated(labels)])
    ev_abort("Duplicate contrast name(s): {.val {dup}}.",
             class = "ev_dup_contrast")
  }
  parts <- Map(function(f, lab) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    df$contrast <- lab
    df
  }, files, labels)
  combined <- dplyr::bind_rows(parts)
  if (!is.null(padj) && padj %in% colnames(combined)) {
    names(combined)[names(combined) == padj] <- "adj.P.Val"
  }
  validated <- ev_validate(combined, x = x, y = y,
                           uniprot = uniprot, species = species)
  idmap_rep <- attr(validated, "ev_idmap_report")
  if (!"adj.P.Val" %in% colnames(validated)) {
    validated <- adjust_p(validated, method = p_adjust)
  }
  if (!is.null(idmap_rep)) {
    attr(validated, "ev_idmap_report") <- idmap_rep
  }
  validated
}
