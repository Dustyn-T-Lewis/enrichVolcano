# Build inst/extdata/idmap/uniprot2symbol_<species>.tsv.gz from the UniProt
# REST stream. Not shipped (data-raw is in .Rbuildignore). Re-runnable to
# refresh the bundled tables.
#
# Output columns: accession (chr), symbol (chr), is_secondary (lgl).
# Primary accessions map to gene_primary (is_secondary = FALSE); each secondary
# accession in sec_acc maps to the same gene_primary (is_secondary = TRUE).

suppressMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
})

dest <- "inst/extdata/idmap"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
dir.create("data-raw/raw", recursive = TRUE, showWarnings = FALSE)

species_tax <- c(human = 9606L, mouse = 10090L, rat = 10116L)

fetch_idmap <- function(species, tax) {
  raw <- file.path("data-raw/raw", paste0("uniprot_idmap_", species, ".tsv"))
  if (!file.exists(raw)) {
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/stream?format=tsv",
      "&query=reviewed:true+AND+organism_id:", tax,
      "&fields=accession,gene_primary,sec_acc"
    )
    message("Downloading UniProt id-map for ", species, " ...")
    utils::download.file(url, raw, quiet = TRUE)
  }
  readr::read_tsv(raw, show_col_types = FALSE)
}

build_one <- function(species, tax) {
  up <- fetch_idmap(species, tax)
  names(up) <- c("accession", "symbol", "sec_acc")
  up <- up |> filter(!is.na(symbol), nzchar(symbol))
  primary <- up |>
    transmute(accession, symbol, is_secondary = FALSE)
  secondary <- up |>
    filter(!is.na(sec_acc), nzchar(sec_acc)) |>
    transmute(accession = sec_acc, symbol) |>
    separate_rows(accession, sep = ";\\s*") |>
    filter(nzchar(accession)) |>
    transmute(accession, symbol, is_secondary = TRUE)
  tab <- bind_rows(primary, secondary) |>
    distinct(accession, .keep_all = TRUE)   # primary wins on collision
  out <- file.path(dest, paste0("uniprot2symbol_", species, ".tsv.gz"))
  readr::write_tsv(tab, out)
  message("Wrote ", out, " (", nrow(tab), " accessions)")
}

for (sp in names(species_tax)) build_one(sp, species_tax[[sp]])
