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
#' one set of column names, so the plots and [run_enrichment()] can read it.
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
#' says which. The signed p-value ranking is the one Reimand et al. (2019)
#' recommend for GSEA.
#'
#' @param x A results table with a contrast column (`contrast` or MSstats'
#'   `Label`), a single-contrast table with `contrast` set, or a named list of
#'   tables whose names are the contrasts.
#' @param contrast Label for a table that has no contrast column.
#' @param species For UniProt accessions, the species whose annotation
#'   package supplies current gene symbols (`"Homo sapiens"`, `"Mus musculus"`,
#'   `"Rattus norvegicus"`). `NULL` keeps the table's own symbols. Tables keyed
#'   by symbols rather than accessions keep theirs either way.
#' @param protein,gene,logfc,t,p,padj,abundance Column names, when detection
#'   does not find them. `logfc` names the fold-change column.
#'
#' @section Gene symbols:
#' Gene sets list gene symbols, so each UniProt accession is mapped to its
#' current symbol through the species' annotation package (isoform suffixes
#' such as `-2` are dropped first). A symbol from the search engine's FASTA can
#' be out of date: `O00483` is `COXFA4`, formerly `NDUFA4`.
#'
#' @return A data frame of class `enrichVolcano_da` with columns `protein`,
#'   `gene`, `contrast`, `logFC`, `t`, `p`, `padj`, `abundance`, `rank` and
#'   `rank_stat`, followed by any input columns it did not use.
#' @references
#' Reimand J, Isserlin R, Voisin V, et al. (2019). Pathway enrichment analysis
#' and visualization of omics data using g:Profiler, GSEA, Cytoscape and
#' EnrichmentMap. Nature Protocols 14:482-517. \doi{10.1038/s41596-018-0103-9}
#' @export
#' @examples
#' tbl <- data.frame(
#'   logFC = c(1.2, -0.8), t = c(4.1, -3.2), P.Value = c(1e-4, 2e-3),
#'   adj.P.Val = c(1e-3, 0.01), row.names = c("P31040", "Q9UBK2")
#' )
#' as_da(tbl, contrast = "Aging")
as_da <- function(x, contrast = NULL, species = "Homo sapiens", protein = NULL, gene = NULL,
                  logfc = NULL, t = NULL, p = NULL, padj = NULL, abundance = NULL) {
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
      refuse_shape(x)
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
  out <- lookup_genes(out, species)
  core <- c("protein", "gene", "contrast", "logFC", "t", "p", "padj", "abundance", "rank", "rank_stat")
  structure(out[c(core, setdiff(names(out), core))], class = c("enrichVolcano_da", "data.frame"))
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
    refuse_shape(tbl)
    ev_abort_missing_column(tbl, "logFC", "logFC", "x")
  }
  stat <- vapply(c("t", "p", "padj", "abundance", "gene", "protein"), field, character(1))
  if (all(is.na(stat[c("t", "p", "padj")]))) {
    ev_abort("Results need one of t, p or padj to rank proteins.", class = "enrichVolcano_column_error")
  }
  column <- function(name) if (is.na(stat[[name]])) NA_real_ else as.numeric(tbl[[stat[[name]]]])

  ids <- protein_ids(tbl, stat[["protein"]], stat[["gene"]])
  res <- data.frame(
    protein = ids$ids,
    logFC = as.numeric(tbl[[logfc_col]]),
    t = column("t"), p = column("p"), padj = column("padj"), abundance = column("abundance"),
    stringsAsFactors = FALSE
  )
  res$gene <- if (is.na(stat[["gene"]])) NA_character_ else unwrap_excel(tbl[[stat[["gene"]]]])
  res$gene[!nzchar(res$gene)] <- NA
  signed <- function(pv) sign(res$logFC) * -log10(pmax(pv, .Machine$double.xmin))
  has_t <- !is.na(stat[["t"]])
  has_p <- !is.na(stat[["p"]])
  res$rank <- if (has_t) res$t else if (has_p) signed(res$p) else signed(res$padj)
  res$rank_stat <- if (has_t) "t" else if (has_p) "signed -log10(p)" else "signed -log10(padj)"
  used <- c(logfc_col, stat[!is.na(stat)], ids$col)
  extras <- tbl[setdiff(names(tbl), used)]
  rownames(res) <- NULL
  rownames(extras) <- NULL
  cbind(res, extras)
}

protein_ids <- function(tbl, col, gene_col) {
  if (!is.na(col)) {
    return(list(ids = as.character(tbl[[col]]), col = col))
  }
  unnamed <- intersect(c("X", "", "...1"), names(tbl))
  if (length(unnamed) > 0) {
    return(list(ids = as.character(tbl[[unnamed[1]]]), col = unnamed[1]))
  }
  if (!all(grepl("^[0-9]+$", rownames(tbl)))) {
    return(list(ids = rownames(tbl), col = NULL))
  }
  if (!is.na(gene_col)) {
    return(list(ids = unwrap_excel(tbl[[gene_col]]), col = gene_col))
  }
  ev_abort_missing_column(tbl, "protein", "protein", "x")
}

lookup_genes <- function(d, species) {
  lead <- sub(";.*", "", d$protein)
  accession <- grepl(uniprot_pattern, lead)
  if (!is.null(species) && any(accession)) {
    symbols <- map_symbols(lead[accession], species)
    n_ids <- length(unique(d$protein[accession]))
    n_mapped <- length(unique(d$protein[accession][!is.na(symbols)]))
    d$gene[accession] <- ifelse(is.na(symbols), d$gene[accession], symbols)
    ev_inform("{n_mapped} of {n_ids} accession{?s} mapped to {species} symbols.",
      class = "enrichVolcano_symbol_lookup"
    )
  }
  fill <- is.na(d$gene) & !accession
  d$gene[fill] <- d$protein[fill]
  d
}

# F-test and wide tables have no single fold change per protein and contrast.
refuse_shape <- function(tbl) {
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
