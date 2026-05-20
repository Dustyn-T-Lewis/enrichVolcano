#' Run enrichment (fgsea + ORA) against one or many databases
#'
#' @details
#' Each entry in `databases` (or each registered name) produces fgsea and/or
#' ORA results. fgsea ranks genes by `rank_by` (default `pi_eq2`). ORA uses
#' significant genes (`pi_eq2 < 0.05` if present, else `adj.P.Val < 0.05`)
#' against the supplied or inferred background.
#'
#' @param data Tidy long tibble (post `pi_score()` / `adjust_p()`).
#' @param contrast Character; one or many contrast names.
#' @param databases Either character vector of registered DB names (resolved
#'   via the registry), a named list of pathway lists, or path(s) to GMT files.
#' @param species Character (for registered DBs).
#' @param enrich_mode `c("fgsea", "ora")` — either or both.
#' @param rank_by Column to rank genes by for fgsea.
#' @param min_size,max_size Gene-set size filters.
#' @param background Optional character vector of background genes for ORA;
#'   default is all genes in `data` for the given contrast.
#' @param exclude_terms Regex of pathway names to exclude (default removes
#'   DISEASE/CANCER/TUMOR terms).
#' @param nperm Permutations for fgsea (default 10000).
#' @param seed Integer seed.
#' @return Tidy tibble: contrast | database | pathway | pval | padj | NES |
#'   size | leading_edge | mode | direction. fgsea-mode rows also carry
#'   `ES` and `log2err` (NA for ORA rows). The pathway list used for
#'   enrichment is attached as `attr(., "ev_pathways")` and the per-contrast
#'   gene-level ranking vectors as `attr(., "ev_stats")`, both required by
#'   `ev_collapse(method = "collapse_fgsea")`.
#' @export
ev_enrich <- function(data, contrast,
                      databases,
                      species = "human",
                      enrich_mode = c("fgsea", "ora"),
                      rank_by = "pi_eq2",
                      min_size = 10, max_size = 500,
                      background = NULL,
                      exclude_terms = "DISEASE|CANCER|TUMOR",
                      nperm = 10000, seed = 42) {
  enrich_mode <- match.arg(enrich_mode, choices = c("fgsea", "ora"),
                           several.ok = TRUE)
  pathway_list <- ev_resolve_databases(databases, species = species)
  results <- list()
  stats_by_contrast <- list()
  for (ctr in contrast) {
    sub <- data[data$contrast == ctr, , drop = FALSE]
    if (nrow(sub) == 0) {
      ev_warn("No rows for contrast {.val {ctr}}; skipping.",
              class = "ev_empty_contrast")
      next
    }
    stats_by_contrast[[ctr]] <- ev_rank_stats(sub, rank_by)
    for (db_name in names(pathway_list)) {
      paths <- pathway_list[[db_name]]
      paths <- paths[!grepl(exclude_terms, names(paths), ignore.case = TRUE)]
      if ("fgsea" %in% enrich_mode) {
        results[[paste(ctr, db_name, "fgsea", sep = "::")]] <-
          ev_fgsea_one(sub, paths, db_name, ctr, rank_by,
                       min_size, max_size, nperm, seed)
      }
      if ("ora" %in% enrich_mode) {
        results[[paste(ctr, db_name, "ora", sep = "::")]] <-
          ev_ora_one(sub, paths, db_name, ctr, background,
                     min_size, max_size)
      }
    }
  }
  out <- dplyr::bind_rows(results)
  if (nrow(out) == 0) {
    out <- tibble::tibble(
      contrast     = character(0),
      database     = character(0),
      pathway      = character(0),
      pval         = double(0),
      padj         = double(0),
      log2err      = double(0),
      ES           = double(0),
      NES          = double(0),
      size         = integer(0),
      leading_edge = character(0),
      mode         = character(0),
      direction    = character(0)
    )
  }
  attr(out, "ev_pathways") <- pathway_list
  attr(out, "ev_stats")    <- stats_by_contrast
  out
}

# Build the gene-level ranking vector fgsea consumes. fgseaMultilevel and
# collapsePathways both require unique gene names; when a gene appears more
# than once (protein isoforms, or test data with replace = TRUE), keep the
# row with the largest absolute rank value — the most informative signal.
ev_rank_stats <- function(sub, rank_by) {
  stats <- stats::setNames(sub[[rank_by]], sub$gene)
  stats <- stats[!is.na(stats)]
  if (anyDuplicated(names(stats))) {
    stats <- tapply(stats, names(stats), function(x) x[which.max(abs(x))])
    stats <- unlist(stats, use.names = TRUE)
  }
  sort(stats, decreasing = TRUE)
}

ev_fgsea_one <- function(sub, paths, db_name, ctr, rank_by,
                          min_size, max_size, nperm, seed) {
  stats <- ev_rank_stats(sub, rank_by)
  withr::with_seed(seed, {
    fg <- fgsea::fgseaMultilevel(
      pathways = paths,
      stats = stats,
      minSize = min_size,
      maxSize = max_size,
      eps = 0
    )
  })
  if (nrow(fg) == 0) return(NULL)
  tibble::tibble(
    contrast = ctr,
    database = db_name,
    pathway = fg$pathway,
    pval = fg$pval,
    padj = fg$padj,
    log2err = fg$log2err,
    ES = fg$ES,
    NES = fg$NES,
    size = fg$size,
    leading_edge = vapply(fg$leadingEdge, paste, character(1),
                           collapse = ";"),
    mode = "fgsea",
    direction = ifelse(fg$NES > 0, "up", "down")
  )
}

ev_ora_one <- function(sub, paths, db_name, ctr, background,
                       min_size, max_size) {
  sig_col <- if ("pi_eq2" %in% colnames(sub)) "pi_eq2" else "adj.P.Val"
  sig <- sub$gene[!is.na(sub[[sig_col]]) & sub[[sig_col]] < 0.05]
  bg <- background %||% sub$gene
  ora <- fgsea::fora(pathways = paths, genes = sig, universe = bg,
                     minSize = min_size, maxSize = max_size)
  if (nrow(ora) == 0) return(NULL)
  tibble::tibble(
    contrast = ctr,
    database = db_name,
    pathway = ora$pathway,
    pval = ora$pval,
    padj = ora$padj,
    NES = NA_real_,
    size = ora$size,
    leading_edge = vapply(ora$overlapGenes, paste, character(1),
                           collapse = ";"),
    mode = "ora",
    direction = NA_character_
  )
}

ev_resolve_databases <- function(databases, species) {
  if (is.list(databases) && !is.character(databases)) return(databases)
  reg <- list_databases()
  resolved <- list()
  for (db_name in databases) {
    row <- reg[reg$name == db_name, ]
    if (nrow(row) == 0) {
      ev_abort("Unknown database {.val {db_name}}.",
               class = "ev_unknown_database")
    }
    species_ok <- grepl(species, row$species, fixed = TRUE)
    if (!species_ok) {
      ev_warn("Database {.val {db_name}} does not list {.val {species}}; skipping.",
              class = "ev_species_unsupported")
      next
    }
    resolved[[db_name]] <- ev_fetch_pathways(db_name, species, row$source)
  }
  resolved
}

ev_fetch_pathways <- function(db_name, species, source) {
  if (grepl("msigdbr", source)) {
    cat_info <- ev_msigdbr_category(db_name)
    df <- msigdbr::msigdbr(
      species = ev_msigdbr_species(species),
      collection = cat_info$collection,
      subcollection = cat_info$subcollection
    )
    split(df$gene_symbol, df$gs_name)
  } else if (source == "bundled_gmt") {
    fn <- paste0(db_name, "_", species, ".gmt")
    path <- system.file("extdata", "gmt", fn, package = "enrichVolcano")
    if (!nzchar(path)) {
      ev_abort("Missing bundled GMT: {.file {fn}}.",
               class = "ev_missing_gmt")
    }
    fgsea::gmtPathways(path)
  } else {
    ev_abort("Fetch-on-demand source {.val {source}} not yet supported.",
             class = "ev_fetch_pending")
  }
}

ev_msigdbr_category <- function(db_name) {
  map <- list(
    hallmark           = list(collection = "H",  subcollection = NULL),
    reactome           = list(collection = "C2", subcollection = "CP:REACTOME"),
    kegg_legacy        = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
    kegg_medicus       = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
    wikipathways       = list(collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
    biocarta           = list(collection = "C2", subcollection = "CP:BIOCARTA"),
    pid                = list(collection = "C2", subcollection = "CP:PID"),
    go_bp              = list(collection = "C5", subcollection = "GO:BP"),
    go_cc              = list(collection = "C5", subcollection = "GO:CC"),
    go_mf              = list(collection = "C5", subcollection = "GO:MF"),
    hpo                = list(collection = "C5", subcollection = "HPO"),
    hallmark_oncogenic = list(collection = "C6", subcollection = NULL),
    immunesigdb        = list(collection = "C7", subcollection = "IMMUNESIGDB"),
    cell_type_sig      = list(collection = "C8", subcollection = NULL),
    mirdb              = list(collection = "C3", subcollection = "MIR:MIRDB"),
    tft_gtrd           = list(collection = "C3", subcollection = "TFT:GTRD")
  )
  out <- map[[db_name]]
  if (is.null(out)) {
    ev_abort("No msigdbr map for {.val {db_name}}.",
             class = "ev_msigdbr_unmapped")
  }
  out
}

ev_msigdbr_species <- function(species) {
  map <- list(
    human = "Homo sapiens", mouse = "Mus musculus",
    rat = "Rattus norvegicus", zebrafish = "Danio rerio",
    fly = "Drosophila melanogaster",
    yeast = "Saccharomyces cerevisiae", pig = "Sus scrofa"
  )
  out <- map[[species]]
  if (is.null(out)) {
    ev_abort("Unknown species {.val {species}}.",
             class = "ev_unknown_species")
  }
  out
}
