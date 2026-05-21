# Script to (re)build bundled GMT files. Not shipped to users.
# Run interactively after manual downloads to ensure license attribution.
# Sources required (manually download into data-raw/raw/):
#   - MitoCarta3 human/mouse/rat .xls files from Broad
#   - HPA subcellular_location.tsv from proteinatlas.org
#   - UniProt subcellular .tsv (one per species) via REST or pre-download

library(readr)
library(dplyr)
library(tidyr)
library(purrr)

dest <- "inst/extdata/gmt"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

# ---- MitoCarta 3.0 ----
build_mitocarta_gmt <- function(input_path, species) {
  if (!file.exists(input_path)) {
    message("Skipping MitoCarta ", species, ": ", input_path, " not found.")
    return(invisible())
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    message("readxl not installed; skipping MitoCarta build.")
    return(invisible())
  }
  paths <- readxl::read_excel(input_path, sheet = "C MitoPathways")
  paths <- paths[!is.na(paths$MitoPathway) & !is.na(paths$Genes), ]
  out <- file.path(dest, paste0("mitocarta3_", species, ".gmt"))
  lines <- vapply(seq_len(nrow(paths)), function(i) {
    genes <- strsplit(paths$Genes[i], ",\\s*")[[1]]
    genes <- trimws(genes)
    genes <- genes[nzchar(genes)]
    name <- paste0("MITOCARTA_", gsub("[^A-Za-z0-9]+", "_", paths$MitoPathway[i]))
    paste(c(name, "MitoCarta3.0", genes), collapse = "\t")
  }, character(1))
  writeLines(lines, out)
  message("Wrote ", out, " (", nrow(paths), " MitoPathways)")
}
build_mitocarta_gmt("data-raw/raw/Human.MitoCarta3.0.xls", "human")
build_mitocarta_gmt("data-raw/raw/Mouse.MitoCarta3.0.xls", "mouse")

# ---- HPA subcellular ----
build_hpa_gmt <- function(input_path) {
  if (!file.exists(input_path)) {
    message("Skipping HPA: ", input_path, " not found.")
    return(invisible())
  }
  hpa <- read_tsv(input_path, show_col_types = FALSE) |>
    select(`Gene name`, `Main location`) |>
    filter(!is.na(`Main location`)) |>
    separate_rows(`Main location`, sep = ";") |>
    group_by(`Main location`) |>
    summarise(genes = list(unique(`Gene name`)), .groups = "drop")
  out <- file.path(dest, "hpa_subcellular_human.gmt")
  lines <- map_chr(seq_len(nrow(hpa)), function(i) {
    paste(c(paste0("HPA_", gsub(" ", "_", hpa$`Main location`[i])),
            "HPA Subcellular", hpa$genes[[i]]), collapse = "\t")
  })
  writeLines(lines, out)
  message("Wrote ", out)
}
build_hpa_gmt("data-raw/raw/subcellular_location.tsv")

# ---- UniProt subcellular ----
# UniProt's REST TSV ships the location as the free-text "Subcellular location
# [CC]" field, e.g.
#   "SUBCELLULAR LOCATION: Cytoplasm {ECO:...}; Multi-pass membrane protein. Nucleus."
# Reduce each entry to its top-level compartment terms: strip the prefix and
# {ECO} evidence tags, drop trailing Note=... text and [Isoform]: markers, split
# into location statements on ".", drop membrane-topology descriptors after the
# first ";", and take the leading compartment before the first ",".
ev_parse_uniprot_cc <- function(x) {
  x <- sub("^SUBCELLULAR LOCATION:\\s*", "", x)
  x <- gsub("\\{[^}]*\\}", "", x)        # ECO evidence tags
  x <- sub("Note=.*$", "", x)            # trailing free-text notes
  x <- gsub("\\[Isoform[^]]*\\]:", "", x) # isoform markers
  stmts <- strsplit(x, "\\.\\s*")[[1]]
  comp <- vapply(stmts, function(s) {
    s <- strsplit(s, ";", fixed = TRUE)[[1]][1]  # drop membrane topology
    s <- strsplit(s, ",", fixed = TRUE)[[1]][1]  # top-level compartment
    trimws(s)
  }, character(1), USE.NAMES = FALSE)
  comp <- comp[!is.na(comp) & nzchar(comp)]
  unique(comp)
}

build_uniprot_gmt <- function(input_path, species) {
  if (!file.exists(input_path)) {
    message("Skipping UniProt ", species, ": ", input_path, " not found.")
    return(invisible())
  }
  up <- read_tsv(input_path, show_col_types = FALSE)
  # Map UniProt REST field names (new = old) to the canonical columns; tolerate
  # files that already use the canonical names.
  rest_map <- c(gene_name = "Gene Names (primary)",
                subcellular_location = "Subcellular location [CC]")
  hit <- rest_map %in% names(up)
  if (any(hit)) up <- dplyr::rename(up, !!!rest_map[hit])
  if (!all(c("gene_name", "subcellular_location") %in% names(up))) {
    stop("UniProt ", species, ": expected columns gene_name + ",
         "subcellular_location (or the UniProt REST equivalents).")
  }
  up <- up |>
    select(gene_name, subcellular_location) |>
    filter(!is.na(gene_name), !is.na(subcellular_location),
           nzchar(gene_name)) |>
    mutate(loc = map(subcellular_location, ev_parse_uniprot_cc)) |>
    select(gene_name, loc) |>
    unnest(loc) |>
    distinct() |>
    group_by(loc) |>
    summarise(genes = list(unique(gene_name)), .groups = "drop")
  out <- file.path(dest, paste0("uniprot_subcellular_", species, ".gmt"))
  lines <- map_chr(seq_len(nrow(up)), function(i) {
    paste(c(paste0("UNIPROT_", gsub(" ", "_", up$loc[i])),
            "UniProt Subcellular", up$genes[[i]]), collapse = "\t")
  })
  writeLines(lines, out)
  message("Wrote ", out, " (", nrow(up), " compartments)")
}
build_uniprot_gmt("data-raw/raw/uniprot_human.tsv", "human")
build_uniprot_gmt("data-raw/raw/uniprot_mouse.tsv", "mouse")
build_uniprot_gmt("data-raw/raw/uniprot_rat.tsv", "rat")

# ---- CORUM 5.0 ----
# Source: data-raw/raw/allComplexes.txt (download from
# https://mips.helmholtz-muenchen.de/corum/download)
build_corum_gmt <- function(input_path = "data-raw/raw/allComplexes.txt") {
  if (!file.exists(input_path)) {
    message("Skipping CORUM: ", input_path, " not found.")
    return(invisible())
  }
  cor <- data.table::fread(input_path, sep = "\t", header = TRUE,
                           fill = TRUE, quote = "\"")
  for (sp_name in c("Human", "Mouse", "Rat")) {
    sub <- cor[cor$organism == sp_name &
               !is.na(cor$subunits_gene_name) &
               cor$subunits_gene_name != "", ]
    if (nrow(sub) == 0) {
      message("No CORUM complexes for ", sp_name)
      next
    }
    out_species <- tolower(sp_name)
    out <- file.path(dest, paste0("corum_", out_species, ".gmt"))
    lines <- vapply(seq_len(nrow(sub)), function(i) {
      genes <- strsplit(sub$subunits_gene_name[i], ";", fixed = TRUE)[[1]]
      genes <- trimws(genes)
      genes <- genes[nzchar(genes)]
      name <- paste0("CORUM_", gsub("[^A-Za-z0-9]+", "_", sub$complex_name[i]))
      paste(c(name, "CORUM_5.0", genes), collapse = "\t")
    }, character(1))
    writeLines(lines, out)
    message("Wrote ", out, " (", nrow(sub), " complexes)")
  }
}
build_corum_gmt()
