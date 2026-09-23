da_candidates <- list(
  protein = c("protein", "Protein", "uniprot_id", "UniProt ID", "feature", "name", "protein_Id", "row"),
  gene = c("gene", "Genes", "gene_symbol", "Gene", "PG.Genes", "Gene.names"),
  contrast = c("contrast", "Label"),
  logFC = c("logFC", "log2FC", "diff", "LFC"),
  t = c("t", "Tvalue", "t_statistic", "statistic"),
  p = c("p", "P.Value", "pvalue", "pval", "p.value"),
  padj = c("padj", "adj.P.Val", "adj.pvalue", "adj_pval", "adjPval", "FDR", "bh"),
  abundance = c("abundance", "AveExpr", "average_intensity", "avg_abundance", "avgAbd")
)

#' Standardise differential-abundance results
#'
#' Reads a DA results table from any common proteomics tool and returns it with
#' one set of column names, so the plots and `run_enrichment()` can read it.
#'
#' @section Recognised tools:
#' limma and limpa `topTable()` (IDs in the row names, or an unnamed first
#' column once written to CSV), proteoDA results and per-contrast CSVs, MSstats
#' `groupComparison()`, proDA `test_diff()`, msqrob2 `topFeatures()`, prolfQua
#' contrast tables and ProtRank output. Any other table works through the
#' column arguments. F-test tables (no fold change) and wide tables with one
#' column per contrast are refused; reshape them to one row per protein and
#' contrast first.
#'
#' @section Ranking statistic:
#' `rank` is the moderated t when the table has one, else
#' `sign(logFC) * -log10(p)`, else `sign(logFC) * -log10(padj)`; `rank_stat`
#' says which.
#'
#' @param x A results table with a contrast column (`contrast` or MSstats'
#'   `Label`), a single-contrast table with `contrast` set, or a named list of
#'   tables whose names are the contrasts.
#' @param contrast Label for a table that has no contrast column.
#' @param protein,gene,logfc,t,p,padj,abundance Column names, when detection
#'   does not find them. `logfc` names the fold-change column.
#'
#' @return A data frame of class `enrichVolcano_da` with columns `protein`,
#'   `gene`, `contrast`, `logFC`, `t`, `p`, `padj`, `abundance`, `rank` and
#'   `rank_stat`.
#' @export
#' @examples
#' tbl <- data.frame(
#'   logFC = c(1.2, -0.8), t = c(4.1, -3.2), P.Value = c(1e-4, 2e-3),
#'   adj.P.Val = c(1e-3, 0.01), row.names = c("P31040", "Q9UBK2")
#' )
#' as_da(tbl, contrast = "Aging")
as_da <- function(x, contrast = NULL, protein = NULL, gene = NULL, logfc = NULL,
                  t = NULL, p = NULL, padj = NULL, abundance = NULL) {
  cols <- list(
    protein = protein, gene = gene, logFC = logfc, t = t, p = p,
    padj = padj, abundance = abundance
  )
  if (is.data.frame(x)) {
    x <- as.data.frame(x)
    found <- intersect(da_candidates$contrast, names(x))
    if (length(found) > 0) {
      names(x)[names(x) == found[1]] <- "contrast"
    } else if (!is.null(contrast)) {
      x$contrast <- contrast
    } else {
      ev_abort(
        c(
          "This table has no contrast column.",
          i = "Set {.arg contrast} to label it, or pass a named list of tables."
        ),
        class = "enrichVolcano_input_error"
      )
    }
  }
  tables <- split_contrasts(x)
  out <- bind_rows_fill(Map(function(tbl, label) {
    tbl$contrast <- NULL
    res <- standardise_da(tbl, cols)
    res$contrast <- label
    res
  }, tables, names(tables)))
  validate_da(out)
  structure(
    out[c("protein", "gene", "contrast", "logFC", "t", "p", "padj", "abundance", "rank", "rank_stat")],
    class = c("enrichVolcano_da", "data.frame")
  )
}

standardise_da <- function(tbl, cols) {
  field <- function(name) {
    col <- cols[[name]]
    if (!is.null(col)) {
      require_columns(tbl, col)
      return(col)
    }
    intersect(da_candidates[[name]], names(tbl))[1]
  }
  logfc_col <- field("logFC")
  if (is.na(logfc_col)) {
    if (any(c("F", "f_statistic") %in% names(tbl))) {
      ev_abort("This looks like an F-test table, which has no per-contrast fold change.",
        class = "enrichVolcano_input_error"
      )
    }
    if (any(grepl("^logFC_", names(tbl)))) {
      ev_abort(
        "This looks like a wide table with one column per contrast; reshape it to one row per protein and contrast.",
        class = "enrichVolcano_input_error"
      )
    }
    ev_abort_missing_column(tbl, "logFC", "logFC", "x")
  }
  stat <- vapply(c("t", "p", "padj", "abundance", "gene", "protein"), field, character(1))
  if (all(is.na(stat[c("t", "p", "padj")]))) {
    ev_abort("Results need one of t, p or padj to rank proteins.", class = "enrichVolcano_column_error")
  }
  column <- function(name) if (is.na(stat[[name]])) NA_real_ else as.numeric(tbl[[stat[[name]]]])

  res <- data.frame(
    protein = protein_ids(tbl, stat[["protein"]]),
    logFC = as.numeric(tbl[[logfc_col]]),
    t = column("t"), p = column("p"), padj = column("padj"), abundance = column("abundance"),
    stringsAsFactors = FALSE
  )
  res$gene <- if (is.na(stat[["gene"]])) NA_character_ else unwrap_excel(tbl[[stat[["gene"]]]])
  ranked <- rank_statistic(res, stat)
  res$rank <- ranked$rank
  res$rank_stat <- ranked$label
  rownames(res) <- NULL
  res
}

protein_ids <- function(tbl, col) {
  if (!is.na(col)) {
    return(as.character(tbl[[col]]))
  }
  unnamed <- intersect(c("X", "", "...1"), names(tbl))
  if (length(unnamed) > 0) {
    return(as.character(tbl[[unnamed[1]]]))
  }
  if (!all(grepl("^[0-9]+$", rownames(tbl)))) {
    return(rownames(tbl))
  }
  ev_abort_missing_column(tbl, "protein", "protein", "x")
}

rank_statistic <- function(res, stat) {
  signed <- function(pv) sign(res$logFC) * -log10(pmax(pv, .Machine$double.xmin))
  if (!is.na(stat[["t"]])) {
    return(list(rank = res$t, label = "t"))
  }
  if (!is.na(stat[["p"]])) {
    return(list(rank = signed(res$p), label = "signed -log10(p)"))
  }
  list(rank = signed(res$padj), label = "signed -log10(padj)")
}

unwrap_excel <- function(x) sub('^="(.*)"$', "\\1", as.character(x))

validate_da <- function(d) {
  if (any(!in_unit_interval(d$p)) || any(!in_unit_interval(d$padj))) {
    ev_abort("`p` and `padj` must lie in [0, 1].", class = "enrichVolcano_data_error")
  }
  n_dup <- sum(duplicated(d[c("contrast", "protein")]))
  if (n_dup > 0) {
    ev_abort("{n_dup} duplicate contrast / protein row{?s}.", class = "enrichVolcano_data_error")
  }
}
