#' Validate and canonicalize input from any supported source
#'
#' Detects input format from column patterns, applies the appropriate
#' adapter, and returns a tidy long tibble with canonical columns
#' (gene, contrast, logFC, P.Value, adj.P.Val).
#'
#' @param data Input data frame from any supported source.
#' @param x Column name for log fold change (defaults to auto-detect).
#' @param y Column name for p-value (defaults to auto-detect).
#' @param lab Column name for gene labels (defaults to auto-detect).
#' @param uniprot Column name holding UniProt accessions (defaults to
#'   auto-detect). When supplied or detected, the input is treated as
#'   UniProt-first and `species` is required.
#' @param species One of `"human"`, `"mouse"`, `"rat"`; required only for
#'   UniProt input. Used to map accessions to gene symbols for enrichment.
#' @return Tidy long tibble with attribute `ev_source` indicating detected format.
#' @export
ev_validate <- function(data, x = NULL, y = NULL, lab = NULL,
                        uniprot = NULL, species = NULL) {
  if (!is.data.frame(data)) {
    ev_abort("`data` must be a data.frame or tibble; got {.cls {class(data)[1]}}.",
             class = "ev_input_not_dataframe")
  }
  src <- ev_detect_source(data)
  data <- switch(src,
    DEP_wide       = ev_pivot_dep(data),
    Perseus_wide   = ev_pivot_perseus(data),
    MaxQuant_wide  = ev_pivot_maxquant(data),
    wide_suffix    = ev_pivot_wide_suffix(data),
    data
  )
  data <- tibble::as_tibble(ev_rename_canonical(data, src, x = x, y = y, lab = lab))
  # MSstats can report log10 fold changes; the rest of the package assumes
  # log2. Rescale when logFC was taken from a log10FC column (and no log2FC
  # column was available to prefer).
  if ("log10FC" %in% colnames(data) && !"log2FC" %in% colnames(data) &&
      is.null(x)) {
    ev_inform("Detected MSstats {.field log10FC}; rescaling to log2 via {.code x * log2(10)}.",
              class = "ev_msstats_log10")
    data$logFC <- data$logFC * log2(10)
  }
  if (!"P.Value" %in% colnames(data) && all(c("t", "df") %in% colnames(data))) {
    ev_inform("No {.field P.Value} column; deriving from {.field t} and {.field df}.",
              class = "ev_derive_p")
    data$P.Value <- 2 * stats::pt(-abs(data$t), df = data$df)
  }
  data <- ev_apply_uniprot(data, uniprot = uniprot, species = species)
  ev_assert_min_cols(data)
  attr(data, "ev_source") <- src
  data
}

# ---- detection ----
ev_detect_source <- function(data) {
  cn <- colnames(data)
  if (any(grepl("sca\\.P\\.Value", cn))) return("DEqMS")
  if (all(c("baseMean", "log2FoldChange", "padj") %in% cn)) return("DESeq2")
  if (all(c("logCPM", "PValue") %in% cn)) return("edgeR")
  if (all(c("Label", "adj.pvalue") %in% cn) &&
      any(c("log2FC", "log10FC") %in% cn)) return("MSstats")
  if (all(c("diff", "pval", "adj_pval", "t_statistic") %in% cn)) return("proDA")
  if (any(grepl("_(diff|p\\.val|p\\.adj)$", cn))) return("DEP_wide")
  if (any(grepl("Student's T-test Difference", cn, fixed = TRUE))) return("Perseus_wide")
  if (any(grepl("^Ratio H/L", cn))) return("MaxQuant_wide")
  if (all(c("logFC", "P.Value", "adj.P.Val") %in% cn) && "average_intensity" %in% cn) return("proteoDA")
  if (all(c("logFC", "P.Value", "adj.P.Val", "B") %in% cn)) return("limma")
  if (any(grepl("^(logFC|P\\.Value|adj\\.P\\.Val|pi_score)_", cn))) return("wide_suffix")
  "tidy_long"
}

# ---- canonical rename ----
ev_col_aliases <- list(
  gene      = c("gene", "Gene", "name", "Protein", "ID", "Genes",
                "Protein.IDs", "Majority protein IDs"),
  contrast  = c("contrast", "comparison", "Label", "coef", "term"),
  logFC     = c("logFC", "log2FoldChange", "log2FC", "log10FC", "diff",
                "Difference"),
  P.Value   = c("P.Value", "PValue", "pvalue", "pval", "p.val", "p_value", "sca.P.Value"),
  adj.P.Val = c("adj.P.Val", "padj", "FDR", "FWER", "adj.pvalue", "adj_pval",
                "qval", "p.adj", "sca.adj.pval")
)

ev_rename_canonical <- function(data, src, x = NULL, y = NULL, lab = NULL) {
  cn <- colnames(data)
  rename_map <- list()
  for (canon in names(ev_col_aliases)) {
    user_arg <- switch(canon, logFC = x, P.Value = y, gene = lab, NULL)
    if (!is.null(user_arg) && user_arg %in% cn) {
      rename_map[[canon]] <- user_arg
    } else {
      hit <- intersect(ev_col_aliases[[canon]], cn)
      if (length(hit) > 0) rename_map[[canon]] <- hit[1]
    }
  }
  for (canon in names(rename_map)) {
    src_col <- rename_map[[canon]]
    if (!identical(src_col, canon)) {
      data[[canon]] <- data[[src_col]]
    }
  }
  data
}

ev_assert_min_cols <- function(data) {
  required <- c("gene", "logFC", "P.Value")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    msg <- cli::format_inline(
      "Missing required column(s): {.field {missing_cols}}."
    )
    ev_abort(msg, class = "ev_input_missing_column")
  }
  invisible(data)
}

# ---- adapters ----
ev_pivot_perseus <- function(data) {
  cn <- colnames(data)
  if ("Gene names" %in% cn) names(data)[names(data) == "Gene names"] <- "gene"
  id_cols <- intersect(colnames(data), c("gene", "Gene", "name", "Protein", "ID"))
  long <- data |>
    tidyr::pivot_longer(
      cols = -tidyselect::all_of(id_cols),
      names_to = c(".value", "contrast"),
      names_pattern = "^(Student's T-test Difference|-Log Student's T-test p-value|Student's T-test q-value) (.+)$"
    ) |>
    dplyr::rename(
      logFC = "Student's T-test Difference",
      `neg_log_p` = "-Log Student's T-test p-value",
      adj.P.Val = "Student's T-test q-value"
    ) |>
    dplyr::mutate(P.Value = 10 ^ (-.data$neg_log_p)) |>
    dplyr::select(-"neg_log_p")
  ev_inform("Detected Perseus {.code -Log p-value}; back-transforming via {.code 10^(-x)}.",
            class = "ev_perseus_backtransform")
  long
}

ev_pivot_maxquant <- function(data) {
  cn <- colnames(data)
  id_cols <- intersect(cn, c("gene", "Gene", "name", "Protein", "ID", "Protein.IDs"))
  long <- data |>
    tidyr::pivot_longer(
      cols = -tidyselect::all_of(id_cols),
      names_to = c(".value", "contrast"),
      names_pattern = "^(Ratio H/L normalized|Ratio H/L significance\\([AB]\\)) (.+)$"
    )
  # Find the resulting column names and rename
  out_cols <- colnames(long)
  if ("Ratio H/L normalized" %in% out_cols) {
    long <- dplyr::rename(long, logFC = "Ratio H/L normalized")
  }
  sig_col <- grep("^Ratio H/L significance", out_cols, value = TRUE)
  if (length(sig_col) > 0) {
    long$P.Value <- long[[sig_col[1]]]
    long[[sig_col[1]]] <- NULL
  }
  long$logFC <- log2(long$logFC)
  long$adj.P.Val <- stats::p.adjust(long$P.Value, method = "BH")
  ev_inform("MaxQuant ratios detected; computing BH adjustment.",
            class = "ev_maxquant_bh")
  long
}

ev_pivot_dep <- function(data) {
  cn <- colnames(data)
  id_cols <- intersect(cn, c("gene", "Gene", "name", "Protein", "ID"))
  long <- data |>
    tidyr::pivot_longer(
      cols = -tidyselect::all_of(id_cols),
      names_to = c("contrast", ".value"),
      names_pattern = "^(.+?)_(diff|p\\.val|p\\.adj|significant)$"
    ) |>
    dplyr::rename(
      logFC = "diff",
      P.Value = "p.val",
      adj.P.Val = "p.adj"
    )
  long
}

ev_pivot_wide_suffix <- function(data) {
  id_cols <- intersect(colnames(data), c("gene", "Gene", "name", "Protein", "ID"))
  prefixes <- c("logFC", "P.Value", "adj.P.Val", "pi_score", "pi_eq2",
                "pi_eq1", "sig_pi", "AveExpr", "t", "B")
  prefix_re <- paste0("^(", paste(gsub("\\.", "\\\\.", prefixes), collapse = "|"), ")_(.+)$")
  keep_cols <- c(id_cols, grep(prefix_re, colnames(data), value = TRUE))
  long <- data[, keep_cols, drop = FALSE] |>
    tidyr::pivot_longer(
      cols = -tidyselect::all_of(id_cols),
      names_to = c(".value", "contrast"),
      names_pattern = prefix_re
    )
  long <- long[!is.na(long$contrast), , drop = FALSE]
  long
}

# ---- UniProt resolution ----
ev_detect_uniprot_col <- function(data) {
  # First pass: column name strongly suggests UniProt accessions
  name_hints <- grepl("^(uniprot|accession|uniprot[._-]?id)$",
                      colnames(data), ignore.case = TRUE)
  for (col in colnames(data)[name_hints]) {
    v <- data[[col]]
    if (!is.character(v) && !is.factor(v)) next
    v <- as.character(v)
    v <- v[!is.na(v) & nzchar(v)]
    if (length(v) == 0) next
    if (mean(ev_is_uniprot_accession(v)) > 0.5) return(col)
  }
  # Second pass: column values are overwhelmingly accession-shaped
  for (col in colnames(data)[!name_hints]) {
    v <- data[[col]]
    if (!is.character(v) && !is.factor(v)) next
    v <- as.character(v)
    v <- v[!is.na(v) & nzchar(v)]
    if (length(v) == 0) next
    if (mean(ev_is_uniprot_accession(v)) > 0.8) return(col)
  }
  NULL
}

ev_apply_uniprot <- function(data, uniprot = NULL, species = NULL) {
  acc_col <- uniprot %||% ev_detect_uniprot_col(data)
  if (is.null(acc_col) || !(acc_col %in% colnames(data))) {
    # legacy symbol-only path
    if (!"uniprot" %in% colnames(data)) data$uniprot <- NA_character_
    if (!"symbol" %in% colnames(data)) {
      data$symbol <- if ("gene" %in% colnames(data)) data$gene else NA_character_
    }
    return(data)
  }
  if (is.null(species)) {
    ev_abort(
      c("UniProt accessions detected in {.field {acc_col}} but {.arg species} is NULL.",
        "i" = "Pass {.code species = \"human\"}, {.val mouse}, or {.val rat}."),
      class = "ev_uniprot_no_species"
    )
  }
  data$uniprot <- as.character(data[[acc_col]])
  map <- ev_map_uniprot(data$uniprot, species)
  data$symbol <- map$symbol
  data$gene <- ifelse(map$mapped, map$symbol, data$uniprot)
  rep <- ev_make_idmap_report(map)
  ev_inform(
    "UniProt mapping: {rep$n_mapped}/{rep$n_input} mapped, {rep$n_unmapped} unmapped.",
    class = "ev_idmap_summary"
  )
  attr(data, "ev_idmap_report") <- rep
  data
}
