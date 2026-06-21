#' Load canonical generic GO slim gene sets
#'
#' Parses the shipped `inst/extdata/goslim_generic.obo` (the official GO
#' Consortium generic slim, frozen for reproducibility — see
#' `inst/extdata/goslim_generic.PROVENANCE.md`) to obtain biological_process
#' slim term IDs, expands each term to include all its BP offspring via
#' `GO.db::GOBPOFFSPRING`, then maps genes to symbols via the species
#' `org.db`. Term IDs come exclusively from the OBO file — no inline list.
#'
#' This replaces the previous build-time `data-raw/build_go_slim.R` which
#' shipped pre-built GMTs from a hardcoded 62-term ID vector. The runtime
#' loader is fully canonical: the GO Consortium publishes the generic slim,
#' and `GO.db` / `org.db` are versioned Bioconductor annotation packages.
#'
#' @param species One of `"human"`, `"mouse"`, or `"rat"`.
#' @param min_size Minimum gene-set size (default `10`).
#' @param max_size Maximum gene-set size (default `500`). Broad parents
#'   (`>500` genes) are too non-specific to be informative in enrichment.
#' @return Named list of character vectors; names are prefixed `GOSLIM_`.
#' @export
load_go_slim <- function(species = c("human", "mouse", "rat"),
                         min_size = 10L,
                         max_size = 500L) {
  species <- match.arg(species)
  if (min_size > max_size)
    ev_abort("{.arg min_size} must be <= {.arg max_size}.",
             class = "ev_goslim_size_args")

  orgdb_pkg <- switch(species,
    human = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db",
    rat   = "org.Rn.eg.db"
  )
  for (pkg in c("GO.db", "AnnotationDbi", orgdb_pkg)) {
    if (!requireNamespace(pkg, quietly = TRUE))
      ev_abort("{.pkg {pkg}} is required for the canonical GO slim loader.",
               class = "ev_missing_goslim_dep")
  }

  obo <- system.file("extdata", "goslim_generic.obo",
                     package = "enrichVolcano")
  if (!nzchar(obo))
    ev_abort("Shipped GO slim OBO not found in installed package.",
             class = "ev_goslim_obo_missing")

  bp_slim <- ev_parse_goslim_obo(obo)
  if (nrow(bp_slim) == 0L)
    ev_abort("No biological_process goslim_generic terms parsed from {.file {obo}}.",
             class = "ev_goslim_obo_empty")

  offspring <- as.list(GO.db::GOBPOFFSPRING)

  orgdb  <- getExportedValue(orgdb_pkg, orgdb_pkg)
  go_tbl <- suppressMessages(AnnotationDbi::select(
    orgdb,
    keys    = AnnotationDbi::keys(orgdb, keytype = "GO"),
    keytype = "GO",
    columns = c("SYMBOL", "ONTOLOGY")))
  go_tbl <- go_tbl[!is.na(go_tbl$ONTOLOGY) & go_tbl$ONTOLOGY == "BP", ]
  go_to_sym <- split(go_tbl$SYMBOL, go_tbl$GO)

  sets <- list()
  for (i in seq_len(nrow(bp_slim))) {
    id <- bp_slim$go_id[i]
    nm <- bp_slim$term[i]
    if (is.na(nm) || !nzchar(nm)) next
    terms <- c(id, offspring[[id]])
    genes <- unique(unlist(go_to_sym[intersect(terms, names(go_to_sym))],
                           use.names = FALSE))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    n <- length(genes)
    if (n < min_size || n > max_size) next
    key <- paste0("GOSLIM_",
                  gsub("^_|_$", "",
                       toupper(gsub("[^A-Za-z0-9]+", "_", nm))))
    sets[[key]] <- genes
  }
  sets
}

# Parse an OBO file and return biological_process terms tagged
# `subset: goslim_generic`. All IDs and names come from the file.
ev_parse_goslim_obo <- function(obo_path) {
  lines <- readLines(obo_path, warn = FALSE)
  term_starts <- which(lines == "[Term]")
  if (length(term_starts) == 0L)
    ev_abort("No {.val [Term]} stanzas in {.file {obo_path}}.",
             class = "ev_goslim_obo_malformed")
  term_ends <- c(term_starts[-1L] - 1L, length(lines))

  ids   <- character(0)
  names <- character(0)
  for (k in seq_along(term_starts)) {
    block <- lines[term_starts[k]:term_ends[k]]
    if (any(block == "is_obsolete: true")) next
    if (!any(grepl("^subset: goslim_generic$", block))) next
    ns <- block[grepl("^namespace: ", block)]
    if (length(ns) == 0L || !grepl("biological_process", ns[1L])) next
    id   <- block[grepl("^id: ",   block)]
    name <- block[grepl("^name: ", block)]
    if (length(id) == 0L || length(name) == 0L) next
    ids   <- c(ids,   sub("^id: ",   "", id[1L]))
    names <- c(names, sub("^name: ", "", name[1L]))
  }
  data.frame(go_id = ids, term = names, stringsAsFactors = FALSE)
}
