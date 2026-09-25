#' Select which volcano points to label
#'
#' @param df Single-contrast DA rows with `protein`, `gene` and `logFC`.
#' @param mode One of `"none"`, `"top_per_direction"`, `"by_significance"`,
#'   `"by_genes"`.
#' @param n Cap for the top modes (per direction for `top_per_direction`).
#' @param rank_by `"significance"` (smaller `p_col` first) or `"logfc"`
#'   (larger `|logFC|` first).
#' @param genes For `"by_genes"`: matched against `gene` and `protein`.
#' @param p_col Column that ranks and filters points.
#' @param p_threshold,logfc_threshold Cutoffs defining the candidate pool.
#' @return The rows to label, with a `label_text` column (the gene symbol, or
#'   the accession when the symbol is missing).
#' @keywords internal
#' @noRd
ev_select_labels <- function(df, mode, n, rank_by, genes,
                             p_col, p_threshold, logfc_threshold) {
  df$label_text <- ifelse(is.na(df$gene) | !nzchar(df$gene), df$protein, df$gene)
  if (mode == "by_genes") {
    return(df[df$gene %in% genes | df$protein %in% genes, , drop = FALSE])
  }
  if (mode == "none") {
    return(df[0, , drop = FALSE])
  }

  sig <- df[[p_col]] < p_threshold & abs(df$logFC) >= logfc_threshold
  pool <- df[sig & !is.na(sig), , drop = FALSE]
  ord <- if (rank_by == "logfc") order(-abs(pool$logFC)) else order(pool[[p_col]])
  pool <- pool[ord, , drop = FALSE]
  if (mode == "by_significance") {
    return(utils::head(pool, n))
  }
  # logFC == 0 belongs to neither direction and is intentionally excluded.
  rbind(utils::head(pool[pool$logFC > 0, , drop = FALSE], n), utils::head(pool[pool$logFC < 0, , drop = FALSE], n))
}

#' Clean a pathway name for display
#'
#' Hallmark and GO slim terms take their label from a table written by hand
#' and shipped with the package (`inst/extdata/term_labels.csv`): `"short"`
#' fits two ring lines, `"clean"` spells every word out. The package's earlier
#' cleaning rules informed the short names. Any other name is cleaned plainly:
#' the database prefix and any `>` or `__` levels above the last are dropped,
#' underscores become spaces, words are title-cased and common acronyms such
#' as DNA and mRNA restored.
#'
#' @param name Character vector of term names.
#' @param width Wrap width in characters. A label that already holds a line
#'   break is not wrapped again.
#' @param style `"short"` or `"clean"`. Only table terms differ between the two.
#' @return Character vector of labels.
#' @export
#' @examples
#' clean_label(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_TCA_CYCLE"))
#' clean_label("HALLMARK_OXIDATIVE_PHOSPHORYLATION", style = "clean", width = 40)
clean_label <- function(name, width = 15, style = c("short", "clean")) {
  style <- rlang::arg_match(style)
  out <- plain_label(name)
  tbl <- term_labels()
  hit <- match(name, tbl$term)
  out[!is.na(hit)] <- tbl[[style]][hit[!is.na(hit)]]
  wrap <- !is.na(out) & !grepl("\n", out, fixed = TRUE)
  out[wrap] <- stringr::str_wrap(out[wrap], width = width)
  out
}

term_labels <- function() {
  tbl <- utils::read.csv(system.file("extdata", "term_labels.csv", package = "enrichVolcano"))
  tbl$clean <- gsub("\\n", "\n", tbl$clean, fixed = TRUE)
  tbl$short <- gsub("\\n", "\n", tbl$short, fixed = TRUE)
  tbl
}

plain_label <- function(name) {
  prefix <- "^(HALLMARK|REACTOME|GOSLIM|KEGG_MEDICUS|KEGG_LEGACY|KEGG|GOBP|GOCC|GOMF|WP|CORUM|MITOCARTA)_"
  out <- sub(".*(>|__)", "", sub(prefix, "", name))
  ev_expand_acronyms(stringr::str_to_title(tolower(gsub("_", " ", out))))
}

#' @keywords internal
#' @noRd
ev_expand_acronyms <- function(x) {
  reps <- c(
    "\\bTca\\b" = "TCA", "\\bMapk(\\d?)\\b" = "MAPK\\1",
    "\\bPtk(\\d?)\\b" = "PTK\\1", "\\bRhoa\\b" = "RhoA",
    "\\bRhob\\b" = "RhoB", "\\bRhoc\\b" = "RhoC",
    "\\bGtpase\\b" = "GTPase", "\\bGtp\\b" = "GTP", "\\bGdp\\b" = "GDP",
    "\\bMrna\\b" = "mRNA", "\\bRna\\b" = "RNA", "\\bDna\\b" = "DNA",
    "\\bAtp\\b" = "ATP", "\\bNadh\\b" = "NADH", "\\bFadh\\b" = "FADH",
    "\\bRos\\b" = "ROS", "\\bIfn\\b" = "IFN", "\\bIl-?(\\d+)\\b" = "IL\\1",
    "\\bPi3k\\b" = "PI3K", "\\bAkt\\b" = "AKT", "\\bMtor\\b" = "mTOR",
    "\\bMtorc1\\b" = "MTORC1", "\\bE2f\\b" = "E2F", "\\bMyc\\b" = "MYC",
    "\\bKras\\b" = "KRAS", "\\bTnfa\\b" = "TNFA", "\\bNfkb\\b" = "NF-kB",
    "\\bG2m\\b" = "G2M", "\\bUv\\b" = "UV", "\\bWnt\\b" = "WNT",
    "\\bStat(\\d)\\b" = "STAT\\1", "\\bJak\\b" = "JAK", "\\bTgf\\b" = "TGF",
    "\\bPirnas?\\b" = "piRNAs", "\\bDgc\\b" = "DGC", "\\bMpc\\b" = "MPC",
    "\\bP53\\b" = "p53", "\\bSlc25a\\b" = "SLC25A",
    "\\bRrna\\b" = "rRNA", "\\bTrna\\b" = "tRNA", "\\bSnrna\\b" = "snRNA",
    "\\bEcm\\b" = "ECM", "\\bUch\\b" = "UCH", "\\bHmox(\\d)\\b" = "HMOX\\1",
    "\\bPcp Ce\\b" = "PCP/CE", "\\bPd L1\\b" = "PD-L1",
    "\\bCd(\\d+)\\b" = "CD\\1", "\\bOxphos\\b" = "OXPHOS",
    "\\bComplex ([IiVv]+)\\b" = "Complex \\U\\1",
    "Trna " = "tRNA "
  )
  for (pat in names(reps)) x <- gsub(pat, reps[[pat]], x, perl = TRUE)
  x
}

# Display names for terms: the user's own where given, clean_label() otherwise.
display_labels <- function(terms, labels, width) {
  out <- clean_label(terms, width = width)
  own <- unname(labels[terms])
  hit <- !is.na(own)
  out[hit] <- ifelse(grepl("\n", own[hit], fixed = TRUE), own[hit], stringr::str_wrap(own[hit], width))
  out
}

check_labels <- function(labels) {
  if (is.null(labels)) {
    return(invisible(NULL))
  }
  if (!is.character(labels) || is.null(names(labels)) || any(!nzchar(names(labels)))) {
    ev_abort("{.arg labels} must be a character vector named by term.", class = "enrichVolcano_param_error")
  }
}
