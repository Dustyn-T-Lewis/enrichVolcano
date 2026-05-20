# (Re)build bundled GO Slim generic-BP gene sets. Not shipped to users.
# Requires Bioconductor annotation packages (Suggests):
#   GO.db, org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, AnnotationDbi
# 62 generic GO Slim BP terms (signaling + nervous-system excluded), ported from
# A_YvO_2026/04_Figures/shared/pathway_utils.R::build_goslim_gene_sets().
dest <- "inst/extdata/gmt"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

bp_slim <- c(
  "GO:0000278", "GO:0000910", "GO:0002181", "GO:0002376", "GO:0003012",
  "GO:0003013", "GO:0003014", "GO:0003016", "GO:0005975", "GO:0006091",
  "GO:0006260", "GO:0006281", "GO:0006310", "GO:0006325", "GO:0006351",
  "GO:0006355", "GO:0006399", "GO:0006457", "GO:0006520", "GO:0006629",
  "GO:0006766", "GO:0006886", "GO:0006913", "GO:0006914", "GO:0006954",
  "GO:0007005", "GO:0007010", "GO:0007018", "GO:0007031", "GO:0007059",
  "GO:0007126", "GO:0007155", "GO:0007163", "GO:0007586", "GO:0009100",
  "GO:0012501", "GO:0016071", "GO:0016192", "GO:0023052", "GO:0030154",
  "GO:0030163", "GO:0030198", "GO:0032200", "GO:0034330", "GO:0042060",
  "GO:0042180", "GO:0042254", "GO:0044782", "GO:0048856", "GO:0048870",
  "GO:0050877", "GO:0051604", "GO:0055085", "GO:0055086", "GO:0061024",
  "GO:0065003", "GO:0071941", "GO:0072659", "GO:0098542", "GO:0098754",
  "GO:0140014", "GO:1901135"
)

build_go_slim_gmt <- function(orgdb_pkg, species) {
  if (!requireNamespace(orgdb_pkg, quietly = TRUE) ||
      !requireNamespace("GO.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    message("Annotation packages missing; skipping ", species)
    return(invisible())
  }
  orgdb <- getExportedValue(orgdb_pkg, orgdb_pkg)
  offspring <- as.list(GO.db::GOBPOFFSPRING)
  suppressMessages({
    go_genes <- AnnotationDbi::select(
      orgdb,
      keys = AnnotationDbi::keys(orgdb, keytype = "GO"),
      keytype = "GO", columns = c("SYMBOL", "ONTOLOGY"))
  })
  go_bp <- go_genes[!is.na(go_genes$ONTOLOGY) & go_genes$ONTOLOGY == "BP", ]
  go_to_symbols <- split(go_bp$SYMBOL, go_bp$GO)
  slim_names <- vapply(bp_slim, function(id)
    tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[id]]),
             error = function(e) NA_character_), character(1))

  lines <- character(0)
  for (i in seq_along(bp_slim)) {
    go_id <- bp_slim[i]; nm <- slim_names[i]
    if (is.na(nm)) next
    all_terms <- c(go_id, offspring[[go_id]])
    genes <- unique(unlist(
      go_to_symbols[intersect(all_terms, names(go_to_symbols))],
      use.names = FALSE))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0) next
    set_name <- paste0("GOSLIM_", toupper(gsub("[^A-Za-z0-9]+", "_", nm)))
    lines <- c(lines, paste(c(set_name, "GO Slim (generic BP)", genes),
                            collapse = "\t"))
  }
  out <- file.path(dest, paste0("go_slim_", species, ".gmt"))
  writeLines(lines, out)
  message("Wrote ", out, " (", length(lines), " sets; GO.db ",
          as.character(packageVersion("GO.db")), ")")
}

build_go_slim_gmt("org.Hs.eg.db", "human")
build_go_slim_gmt("org.Mm.eg.db", "mouse")
build_go_slim_gmt("org.Rn.eg.db", "rat")
