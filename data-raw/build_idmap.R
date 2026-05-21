# Build inst/extdata/idmap/uniprot2symbol_<species>.tsv.gz. Not shipped
# (data-raw is in .Rbuildignore). Re-runnable to refresh the bundled tables.
#
# Output columns: accession (chr), symbol (chr), is_secondary (lgl).
# Primary accessions come from the reviewed UniProt REST stream
# (accession + gene_primary, is_secondary = FALSE). Secondary accessions come
# from UniProt's sec_ac.txt index: a secondary accession whose current primary
# is reviewed in this species inherits that primary's symbol (is_secondary =
# TRUE). The REST API has no secondary-accession field, so sec_ac.txt is the
# authoritative reproducible source.

suppressMessages({
  library(readr)
  library(dplyr)
})

if (!file.exists("DESCRIPTION")) {
  stop("Run build_idmap.R from the package root (DESCRIPTION not found).")
}

dest <- "inst/extdata/idmap"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
dir.create("data-raw/raw", recursive = TRUE, showWarnings = FALSE)

species_tax <- c(human = 9606L, mouse = 10090L, rat = 10116L)

download_or_die <- function(url, dest_file) {
  status <- utils::download.file(url, dest_file, quiet = TRUE)
  if (status != 0L || !file.exists(dest_file)) {
    if (file.exists(dest_file)) file.remove(dest_file)
    stop("download.file failed (status ", status, ") for ", url)
  }
}

fetch_primary <- function(species, tax) {
  raw <- file.path("data-raw/raw", paste0("uniprot_idmap_", species, ".tsv"))
  if (!file.exists(raw)) {
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/stream?format=tsv",
      "&query=reviewed:true+AND+organism_id:", tax,
      "&fields=accession,gene_primary"
    )
    message("Downloading UniProt primary id-map for ", species, " ...")
    download_or_die(url, raw)
  }
  up <- readr::read_tsv(raw, show_col_types = FALSE)
  if (ncol(up) != 2L) {
    stop("Expected 2 columns (accession, gene_primary) for ", species,
         "; got ", ncol(up), ". UniProt field layout may have changed.")
  }
  names(up) <- c("accession", "symbol")
  up |> filter(!is.na(symbol), nzchar(symbol))
}

fetch_sec_ac <- function() {
  raw <- "data-raw/raw/sec_ac.txt"
  if (!file.exists(raw)) {
    message("Downloading UniProt sec_ac.txt ...")
    download_or_die(
      "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/sec_ac.txt",
      raw
    )
  }
  lines <- readLines(raw)
  start <- grep("^_+[[:space:]]+_+[[:space:]]*$", lines)
  if (length(start) == 0) stop("sec_ac.txt: could not find the data delimiter.")
  body <- lines[(start[1] + 1):length(lines)]
  body <- body[nzchar(trimws(body))]
  m <- regmatches(body, regexec("^(\\S+)[[:space:]]+(\\S+)[[:space:]]*$", body))
  ok <- lengths(m) == 3L
  data.frame(
    secondary = vapply(m[ok], `[`, character(1), 2L),
    primary   = vapply(m[ok], `[`, character(1), 3L),
    stringsAsFactors = FALSE
  )
}

build_one <- function(species, tax, sec_ac) {
  up <- fetch_primary(species, tax)
  primary <- up |> transmute(accession, symbol, is_secondary = FALSE)
  sec <- sec_ac[sec_ac$primary %in% up$accession, , drop = FALSE]
  sec$symbol <- up$symbol[match(sec$primary, up$accession)]
  secondary <- tibble::tibble(accession = sec$secondary, symbol = sec$symbol,
                              is_secondary = TRUE)
  combined <- bind_rows(primary, secondary)
  tab <- combined |> distinct(accession, .keep_all = TRUE)
  n_ambig <- nrow(combined) - nrow(tab)
  if (n_ambig > 0L) {
    message("  ", species, ": ", n_ambig,
            " duplicate accession rows dropped ",
            "(multi-primary secondaries / primary collisions; first kept)")
  }
  out <- file.path(dest, paste0("uniprot2symbol_", species, ".tsv.gz"))
  readr::write_tsv(tab, out)
  message("Wrote ", out, " (", nrow(tab), " accessions)")
}

sec_ac <- fetch_sec_ac()
for (sp in names(species_tax)) build_one(sp, species_tax[[sp]], sec_ac)
