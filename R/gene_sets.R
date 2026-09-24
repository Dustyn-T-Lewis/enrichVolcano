species_org <- c(
  "Homo sapiens" = "org.Hs.eg.db",
  "Mus musculus" = "org.Mm.eg.db",
  "Rattus norvegicus" = "org.Rn.eg.db"
)

msigdb_collections <- list(
  Hallmark = list(collection = "H", subcollection = NULL),
  Reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  KEGG = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  "GO:BP" = list(collection = "C5", subcollection = "GO:BP")
)

uniprot_pattern <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$"

org_package <- function(species) {
  if (!isTRUE(species %in% names(species_org))) {
    ev_abort(
      "{.arg species} must be one of {.val {names(species_org)}}.",
      class = "enrichVolcano_param_error"
    )
  }
  pkg <- species_org[[species]]
  suppressPackageStartupMessages(
    rlang::check_installed(c("AnnotationDbi", pkg), reason = "to map proteins and gene sets to symbols.")
  )
  pkg
}

map_symbols <- function(ids, species) {
  pkg <- org_package(species)
  db <- getExportedValue(pkg, pkg)
  keys <- sub("-[0-9]+$", "", ids)
  valid <- intersect(unique(keys), AnnotationDbi::keys(db, "UNIPROT"))
  symbols <- if (length(valid) == 0) {
    character(0)
  } else {
    suppressMessages(AnnotationDbi::mapIds(db, valid, "SYMBOL", "UNIPROT", multiVals = "first"))
  }
  stats::setNames(unname(symbols[keys]), ids)
}

#' Load gene-set collections for enrichment
#'
#' Fetches collections at the versions installed on this machine and records
#' those versions, so an analysis can state exactly what it tested.
#'
#' * `"Hallmark"`, `"Reactome"`, `"KEGG"` (KEGG MEDICUS) and `"GO:BP"` come from
#'   msigdbr; mouse and rat sets are human sets mapped through orthologs.
#' * `"GO Slim"` is the biological-process part of the GO Consortium's generic
#'   slim, pinned in this package (release 2026-07-26), with each term's genes
#'   taken from the species' annotation package, counting genes annotated to
#'   the term or any of its descendants. Sets are named like
#'   `GOSLIM_PROTEIN_FOLDING`.
#'
#' Each collection is kept separate so [run_enrichment()] corrects p-values
#' within it.
#'
#' @references
#' Liberzon A, Birger C, Thorvaldsdottir H, et al. (2015). The Molecular
#' Signatures Database hallmark gene set collection. Cell Systems 1(6):417-425.
#' \doi{10.1016/j.cels.2015.12.004}
#'
#' @param databases Any of `"Hallmark"`, `"GO Slim"`, `"Reactome"`, `"KEGG"`,
#'   `"GO:BP"`.
#' @param species `"Homo sapiens"`, `"Mus musculus"` or `"Rattus norvegicus"`.
#' @param min_size,max_size Keep sets with at least `min_size` and at most
#'   `max_size` genes.
#' @return A named list of collections, each a named list of gene symbols, with
#'   a `versions` attribute (msigdbr, GO slim release, annotation package,
#'   species).
#' @export
load_gene_sets <- function(databases = c("Hallmark", "GO Slim"), species = "Homo sapiens",
                           min_size = 15, max_size = 500) {
  known <- c(names(msigdb_collections), "GO Slim")
  unknown <- setdiff(databases, known)
  if (length(unknown) > 0) {
    ev_abort(
      c("Unknown database{?s}: {.val {unknown}}.", i = "Available: {.val {known}}."),
      class = "enrichVolcano_param_error"
    )
  }
  pkg <- org_package(species)
  versions <- list(species = species, annotation = paste(pkg, utils::packageVersion(pkg)))
  sets <- lapply(databases, function(name) {
    if (name == "GO Slim") {
      slim <- go_slim_sets(pkg)
      versions$go_slim <<- attr(slim, "data_version")
      return(slim)
    }
    rlang::check_installed("msigdbr", reason = "to load MSigDB collections.")
    versions$msigdbr <<- as.character(utils::packageVersion("msigdbr"))
    def <- msigdb_collections[[name]]
    tbl <- suppressMessages(msigdbr::msigdbr(
      species = species, collection = def$collection, subcollection = def$subcollection
    ))
    lapply(split(tbl$gene_symbol, tbl$gs_name), unique)
  })
  names(sets) <- databases
  sets <- lapply(sets, function(s) s[lengths(s) >= min_size & lengths(s) <= max_size])
  attr(sets, "versions") <- versions
  sets
}

go_slim_sets <- function(pkg) {
  obo <- readLines(system.file("extdata", "goslim_generic.obo.gz", package = "enrichVolcano"))
  terms <- go_slim_terms(obo)
  db <- getExportedValue(pkg, pkg)
  go2genes <- getExportedValue(pkg, sub("\\.db$", "GO2ALLEGS", pkg))
  terms <- terms[terms$id %in% AnnotationDbi::keys(go2genes), ]
  entrez <- AnnotationDbi::mget(terms$id, go2genes)
  sets <- lapply(entrez, function(genes) {
    symbols <- suppressMessages(AnnotationDbi::mapIds(db, unique(genes), "SYMBOL", "ENTREZID"))
    unique(unname(symbols[!is.na(symbols)]))
  })
  names(sets) <- terms$set_name
  attr(sets, "data_version") <- sub("^data-version: ", "", grep("^data-version:", obo, value = TRUE)[1])
  sets
}

go_slim_terms <- function(obo) {
  stanzas <- split(obo, cumsum(obo == "[Term]"))[-1]
  field <- function(lines, tag) sub(paste0("^", tag, ": "), "", grep(paste0("^", tag, ": "), lines, value = TRUE)[1])
  terms <- data.frame(
    id = vapply(stanzas, field, character(1), tag = "id"),
    name = vapply(stanzas, field, character(1), tag = "name"),
    namespace = vapply(stanzas, field, character(1), tag = "namespace"),
    obsolete = vapply(stanzas, function(lines) any(lines == "is_obsolete: true"), logical(1))
  )
  terms <- terms[terms$namespace == "biological_process" & !terms$obsolete, ]
  terms$set_name <- paste0("GOSLIM_", gsub("^_|_$", "", toupper(gsub("[^A-Za-z0-9]+", "_", terms$name))))
  terms
}
