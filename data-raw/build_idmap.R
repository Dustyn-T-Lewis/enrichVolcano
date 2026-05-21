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
  library(readr); library(dplyr)
})

dest <- "inst/extdata/idmap"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
dir.create("data-raw/raw", recursive = TRUE, showWarnings = FALSE)

species_tax <- c(human = 9606L, mouse = 10090L, rat = 10116L)

fetch_primary <- function(species, tax) {
  raw <- file.path("data-raw/raw", paste0("uniprot_idmap_", species, ".tsv"))
  if (!file.exists(raw)) {
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/stream?format=tsv",
      "&query=reviewed:true+AND+organism_id:", tax,
      "&fields=accession,gene_primary"
    )
    message("Downloading UniProt primary id-map for ", species, " ...")
    utils::download.file(url, raw, quiet = TRUE)
  }
  up <- readr::read_tsv(raw, show_col_types = FALSE)
  names(up) <- c("accession", "symbol")
  up |> filter(!is.na(symbol), nzchar(symbol))
}

fetch_sec_ac <- function() {
  raw <- "data-raw/raw/sec_ac.txt"
  if (!file.exists(raw)) {
    message("Downloading UniProt sec_ac.txt ...")
    utils::download.file(
      "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/sec_ac.txt",
      raw, quiet = TRUE
    )
  }
  lines <- readLines(raw)
  start <- grep("^_+\\s+_+\\s*$", lines)
  if (length(start) == 0) stop("sec_ac.txt: could not find the data delimiter.")
  body <- lines[(start[1] + 1):length(lines)]
  body <- body[nzchar(trimws(body))]
  m <- regmatches(body, regexec("^(\\S+)\\s+(\\S+)\\s*$", body))
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
  tab <- bind_rows(primary, secondary) |>
    distinct(accession, .keep_all = TRUE)
  out <- file.path(dest, paste0("uniprot2symbol_", species, ".tsv.gz"))
  readr::write_tsv(tab, out)
  message("Wrote ", out, " (", nrow(tab), " accessions)")
}

sec_ac <- fetch_sec_ac()
for (sp in names(species_tax)) build_one(sp, species_tax[[sp]], sec_ac)
