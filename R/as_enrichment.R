#' Convert enrichment results into an `enrichment` object
#'
#' The one entry point for results from any enrichment tool. Hand it one
#' table per contrast and it returns a validated [enrichment] that every plot
#' in the package reads.
#'
#' @section Which test produced the table:
#' The test is recognised from the columns it returns:
#' * **fgsea**: `pathway, pval, padj, NES, size, leadingEdge`. Score is NES.
#' * **fry** (`PValue.Mixed`), **mroast** (`PropUp`), **camera** estimating
#'   its own inter-gene correlation (`Correlation`): `NGenes, Direction,
#'   PValue, FDR`, set names in the row names. None of these report an effect
#'   size or a leading edge, so the score is signed \eqn{-\log_{10}}(FDR) and
#'   leading edges are empty.
#' * **camera** with its default fixed correlation and **cameraPR** return
#'   identical columns. Say which with `enrichment_test = "camera"` or
#'   `"cameraPR"`.
#' * **custom**: any table, mapped with the column arguments. Direction comes
#'   from the `direction` column when present, otherwise from the score sign.
#'
#' fry and mroast are self-contained tests (is this set changed at all?);
#' camera, cameraPR and fgsea are competitive (is it changed more than the
#' genes outside it?). The test is recorded so figures can say which.
#'
#' @param x A long table with a `contrast` column, or a named list of tables,
#'   one per contrast, whose names become the contrasts.
#' @param enrichment_test The method that produced `x`: one of `"fgsea"`,
#'   `"fry"`, `"mroast"`, `"camera"`, `"cameraPR"`, `"gseaResult"`, `"ora"`,
#'   `"custom"`. `NULL` recognises it from the columns.
#' @param database A label for every row, such as `"Hallmark"`, when the input
#'   has no `database` column.
#' @param term,score,padj,p,size,direction,leading_edge Column names, used only
#'   when `enrichment_test` is `"custom"` or `"ora"`. `p`, `size` and
#'   `leading_edge` are optional; a leading-edge column may be a list or a
#'   `;`- or `/`-separated string.
#' @param score_type Axis and legend label for a custom score, e.g. `"NES"`.
#'
#' @return An [enrichment] object.
#' @export
#' @examples
#' path <- system.file("extdata", "examples", "yvo_fgsea.csv.gz",
#'   package = "enrichVolcano"
#' )
#' ex <- as_enrichment(read.csv(path))
#' ex
as_enrichment <- function(x, enrichment_test = NULL, database = NULL,
                          term = "term", score = "score", padj = "padj",
                          p = "p", size = "size", direction = "direction",
                          leading_edge = NULL, score_type = NULL) {
  if (!is.null(enrichment_test) && !isTRUE(enrichment_test %in% enrichment_tests)) {
    ev_abort(
      "{.arg enrichment_test} must be one of {.val {enrichment_tests}}.",
      class = "enrichVolcano_param_error"
    )
  }
  tables <- split_contrasts(x)
  test <- unique(vapply(tables, detect_test, character(1), enrichment_test = enrichment_test))
  if (length(test) > 1) {
    ev_abort(
      "Every contrast must come from the same test, not {.val {test}}.",
      class = "enrichVolcano_input_error"
    )
  }
  if (is.null(enrichment_test)) {
    ev_inform("Reading {.val {test}} results.", class = "enrichVolcano_detected_test")
  }
  if (!is.null(score_type) && test != "custom") {
    ev_abort(
      "{.arg score_type} applies only to {.val custom} tables.",
      class = "enrichVolcano_param_error"
    )
  }

  cols <- list(
    term = term, score = score, padj = padj, p = p, size = size,
    direction = direction, leading_edge = leading_edge
  )
  results <- bind_rows_fill(Map(function(tbl, contrast) {
    out <- convert_table(tbl, test, cols)
    out$contrast <- contrast
    stamp_database(out, tbl, database)
  }, tables, names(tables)))
  results <- drop_unscored(results)

  core <- names(results_columns)
  enrichment(
    results = results[c(core, setdiff(names(results), core))],
    metadata = list(
      enrichment_test = test,
      score_type = score_type %||% default_score_type[[test]],
      dedup = NULL
    )
  )
}

default_score_type <- list(
  fgsea = "NES", gseaResult = "NES", custom = "score",
  fry = "signed -log10(FDR)", mroast = "signed -log10(FDR)",
  camera = "signed -log10(FDR)", cameraPR = "signed -log10(FDR)",
  ora = "signed -log10(FDR)"
)

split_contrasts <- function(x) {
  if (is.data.frame(x)) {
    if (!"contrast" %in% names(x)) {
      ev_abort(
        c(
          "Name the contrast these results belong to.",
          i = "Pass a named list, e.g. {.code as_enrichment(list(Aging = res))},
               or a table with a {.field contrast} column."
        ),
        class = "enrichVolcano_input_error"
      )
    }
    x <- as.data.frame(x)
    return(split(x, factor(x$contrast, levels = unique(x$contrast))))
  }
  if (!is.list(x) || is.object(x)) {
    ev_abort(
      "{.arg x} must be a data frame or a named list of them.",
      class = "enrichVolcano_input_error"
    )
  }
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    ev_abort(
      "Every table in {.arg x} needs a contrast name, e.g. {.code list(Aging = res)}.",
      class = "enrichVolcano_input_error"
    )
  }
  lapply(x, function(tbl) {
    if (!is.data.frame(tbl)) {
      ev_abort(
        "Each element of {.arg x} must be a results table, not {.cls {class(tbl)[1]}}.",
        class = "enrichVolcano_input_error"
      )
    }
    as.data.frame(tbl)
  })
}

detect_test <- function(tbl, enrichment_test = NULL) {
  if (!is.null(enrichment_test)) {
    return(enrichment_test)
  }
  cols <- names(tbl)
  if (all(c("NES", "leadingEdge", "pval") %in% cols)) {
    return("fgsea")
  }
  if (all(c("NGenes", "Direction", "FDR") %in% cols)) {
    if ("PropUp" %in% cols) {
      return("mroast")
    }
    if ("PValue.Mixed" %in% cols) {
      return("fry")
    }
    if ("Correlation" %in% cols) {
      return("camera")
    }
    ev_abort(
      c(
        "These columns fit both {.val camera} and {.val cameraPR}.",
        i = "Set {.arg enrichment_test} to say which produced them."
      ),
      class = "enrichVolcano_input_error"
    )
  }
  ev_abort(
    c(
      "Can't tell which test produced this table.",
      i = "Set {.arg enrichment_test}; use {.val custom} with the column arguments
           for anything else."
    ),
    class = "enrichVolcano_input_error"
  )
}

convert_table <- function(tbl, test, cols) {
  switch(test,
    fgsea = from_fgsea(tbl),
    fry = ,
    mroast = ,
    camera = ,
    cameraPR = from_limma(tbl),
    from_table(tbl, test, cols)
  )
}

from_fgsea <- function(tbl) {
  require_columns(tbl, c("pathway", "pval", "padj", "NES", "size", "leadingEdge"))
  new_results(tbl,
    term = tbl$pathway, score = tbl$NES, p = tbl$pval, padj = tbl$padj,
    size = tbl$size, leading_edge = split_genes(tbl$leadingEdge),
    used = c("pathway", "pval", "padj", "NES", "size", "leadingEdge")
  )
}

from_limma <- function(tbl) {
  require_columns(tbl, c("NGenes", "Direction", "PValue", "FDR"))
  direction <- tolower(tbl$Direction)
  new_results(tbl,
    term = rownames(tbl), score = signed_log_fdr(tbl$FDR, direction),
    p = tbl$PValue, padj = tbl$FDR, size = tbl$NGenes,
    direction = direction, used = c("NGenes", "Direction", "PValue", "FDR")
  )
}

from_table <- function(tbl, test, cols) {
  needed <- c(cols$term, cols$padj, if (test == "custom") cols$score, if (test == "ora") cols$direction)
  require_columns(tbl, needed)
  optional <- function(col) if (!is.null(col) && col %in% names(tbl)) tbl[[col]] else NA_real_
  direction <- if (cols$direction %in% names(tbl)) tolower(tbl[[cols$direction]])
  score <- if (test == "ora") signed_log_fdr(tbl[[cols$padj]], direction) else tbl[[cols$score]]
  new_results(tbl,
    term = as.character(tbl[[cols$term]]), score = score,
    p = optional(cols$p), padj = tbl[[cols$padj]], size = optional(cols$size),
    direction = direction,
    leading_edge = if (!is.null(cols$leading_edge)) split_genes(tbl[[cols$leading_edge]]),
    used = unlist(cols[c("term", "score", "padj", "p", "size", "direction", "leading_edge")])
  )
}

new_results <- function(tbl, term, score, p, padj, size, used,
                        direction = NULL, leading_edge = NULL) {
  out <- data.frame(
    term = as.character(term), score = as.numeric(score),
    p = as.numeric(p), padj = as.numeric(padj), size = as.numeric(size),
    direction = direction %||% ifelse(score > 0, "up", "down"),
    stringsAsFactors = FALSE
  )
  out$leading_edge <- leading_edge %||% rep(list(character(0)), nrow(out))
  extras <- tbl[setdiff(names(tbl), c(used, "contrast", "database"))]
  rownames(extras) <- NULL
  cbind(out, extras)
}

signed_log_fdr <- function(fdr, direction) {
  ifelse(direction == "up", 1, -1) * -log10(pmax(fdr, .Machine$double.xmin))
}

split_genes <- function(x) {
  if (is.list(x)) {
    return(lapply(x, as.character))
  }
  lapply(strsplit(as.character(x), "[;/]"), function(g) g[!is.na(g) & nzchar(g)])
}

require_columns <- function(tbl, cols) {
  for (col in cols) {
    if (!col %in% names(tbl)) {
      ev_abort_missing_column(tbl, col, col, "x")
    }
  }
}

stamp_database <- function(out, tbl, database) {
  if ("database" %in% names(tbl)) {
    if (!is.null(database)) {
      ev_abort(
        "{.arg database} labels unlabelled input, but this table already has a {.field database} column.",
        class = "enrichVolcano_input_error"
      )
    }
    out$database <- as.character(tbl$database)
  } else {
    out$database <- database %||% NA_character_
  }
  out
}

drop_unscored <- function(results) {
  unscored <- is.na(results$score)
  if (any(unscored)) {
    ev_inform("Dropped {sum(unscored)} row{?s} with no score.", class = "enrichVolcano_dropped_rows")
  }
  results <- results[!unscored, , drop = FALSE]
  rownames(results) <- NULL
  results
}

bind_rows_fill <- function(tables) {
  cols <- unique(unlist(lapply(tables, names)))
  tables <- lapply(tables, function(tbl) {
    tbl[setdiff(cols, names(tbl))] <- NA
    tbl[cols]
  })
  out <- do.call(rbind, unname(tables))
  rownames(out) <- NULL
  out
}
