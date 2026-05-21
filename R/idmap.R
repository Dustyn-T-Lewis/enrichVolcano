# UniProt accession -> gene symbol mapping, backed by bundled per-species tables
# built in data-raw/build_idmap.R.

.ev_idmap_cache <- new.env(parent = emptyenv())

#' Load the bundled UniProt-to-symbol table for a species
#'
#' @param species One of `"human"`, `"mouse"`, `"rat"`.
#' @return Data frame with columns `accession`, `symbol`, `is_secondary`.
#' @keywords internal
#' @noRd
ev_load_idmap <- function(species) {
  cached <- .ev_idmap_cache[[species]]
  if (!is.null(cached)) return(cached)
  fn <- paste0("uniprot2symbol_", species, ".tsv.gz")
  path <- system.file("extdata", "idmap", fn, package = "enrichVolcano")
  if (!nzchar(path)) {
    ev_abort(
      "No UniProt id-map for {.val {species}} (supported: human, mouse, rat).",
      class = "ev_no_idmap"
    )
  }
  tab <- data.table::fread(
    path, sep = "\t", header = TRUE, data.table = FALSE,
    colClasses = list(character = c("accession", "symbol"),
                      logical = "is_secondary")
  )
  .ev_idmap_cache[[species]] <- tab
  tab
}

#' Map UniProt accessions to gene symbols
#'
#' Strips isoform suffixes (`-\\d+`) for lookup only; the returned `uniprot`
#' column keeps the original accession. Secondary accessions resolve to the
#' primary entry's symbol and are flagged.
#'
#' @param accessions Character vector of UniProt accessions.
#' @param species One of `"human"`, `"mouse"`, `"rat"`.
#' @return Tibble: `uniprot`, `symbol`, `mapped`, `isoform_stripped`,
#'   `via_secondary`.
#' @keywords internal
#' @noRd
ev_map_uniprot <- function(accessions, species) {
  tab <- ev_load_idmap(species)
  orig <- as.character(accessions)
  lookup <- sub("-\\d+$", "", orig)
  idx <- match(lookup, tab$accession)
  tibble::tibble(
    uniprot = orig,
    symbol = tab$symbol[idx],
    mapped = !is.na(idx),
    isoform_stripped = lookup != orig,
    via_secondary = !is.na(idx) & tab$is_secondary[idx]
  )
}

#' Summarise a mapping result
#' @keywords internal
#' @noRd
ev_make_idmap_report <- function(map_tbl) {
  list(
    n_input = nrow(map_tbl),
    n_mapped = sum(map_tbl$mapped),
    n_isoform_stripped = sum(map_tbl$isoform_stripped),
    n_via_secondary = sum(map_tbl$via_secondary, na.rm = TRUE),
    n_unmapped = sum(!map_tbl$mapped),
    unmapped = map_tbl$uniprot[!map_tbl$mapped]
  )
}

#' Retrieve the UniProt mapping report attached to a validated table
#'
#' @param x A tibble produced by `ev_validate()` / `ev_read_contrasts()` from
#'   UniProt input.
#' @return The mapping report list, or `NULL` if the input was not UniProt-based.
#' @export
ev_idmap_report <- function(x) {
  attr(x, "ev_idmap_report")
}

#' Test whether strings look like UniProt accessions
#'
#' Matches the Swiss-Prot accession grammar, allowing an optional `-N` isoform
#' suffix.
#' @keywords internal
#' @noRd
ev_is_uniprot_accession <- function(x) {
  grepl(
    "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$",
    x
  )
}
