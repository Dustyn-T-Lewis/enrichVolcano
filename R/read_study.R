study_sheets <- c("matrix", "samples", "contrasts", "da_results", "weights")

#' Read a study from a workbook or a folder of CSV files
#'
#' A study is up to five tables, linked by name:
#' * `da_results` (required): differential-abundance results in any format
#'   [as_da()] reads, with a contrast column.
#' * `matrix`: one row per protein; the first column holds the protein IDs and
#'   every other column is a sample.
#' * `samples`: `sample` and `group`, optionally `subject` for repeated
#'   measures, plus any covariate columns.
#' * `contrasts`: `name` and `expression`, arithmetic on group names such as
#'   `Post - Pre` or `(B_Post - B_Pre) - (A_Post - A_Pre)`.
#' * `weights` (optional): precision weights shaped like `matrix`.
#'
#' `da_results` alone is enough for fgsea and the figures; fry and camera also
#' need `matrix`, `samples` and `contrasts`. Every link is checked by name:
#' matrix columns against `samples$sample`, contrast expressions against
#' `samples$group`, `da_results` contrasts against `contrasts$name`, and
#' `weights` against `matrix`. Group names must be syntactic (letters, digits,
#' `.` and `_`), because `-` inside a name would read as subtraction.
#'
#' @param path An `.xlsx` workbook with sheets named as above, or a folder of
#'   `<sheet>.csv` or `<sheet>.csv.gz` files.
#' @param species Passed to [as_da()] for the gene-symbol lookup.
#' @return A list of class `enrichVolcano_study` with `da`, `matrix`,
#'   `samples`, `contrasts` and `weights` (absent parts are `NULL`).
#' @export
read_study <- function(path, species = "Homo sapiens") {
  sheets <- read_sheets(path)
  if (is.null(sheets$da_results)) {
    ev_abort("A study needs a {.field da_results} sheet.", class = "enrichVolcano_input_error")
  }
  design_parts <- c("matrix", "samples", "contrasts")
  present <- design_parts[design_parts %in% names(sheets)]
  if (length(present) %in% 1:2) {
    ev_abort(
      "fry and camera need matrix, samples and contrasts together; missing {.field {setdiff(design_parts, present)}}.",
      class = "enrichVolcano_input_error"
    )
  }
  if (!is.null(sheets$weights) && length(present) == 0) {
    ev_abort("A {.field weights} sheet needs matrix, samples and contrasts.", class = "enrichVolcano_input_error")
  }

  study <- list(
    da = as_da(sheets$da_results, species = species),
    matrix = NULL, samples = NULL, contrasts = NULL, weights = NULL
  )
  if (length(present) == 3) {
    study$samples <- check_samples(sheets$samples)
    study$contrasts <- check_contrasts(sheets$contrasts, study$samples$group, study$da$contrast)
    study$matrix <- numeric_matrix(sheets$matrix, "matrix")
    check_sample_links(study$matrix, study$samples)
    study$matrix <- study$matrix[, study$samples$sample, drop = FALSE]
    if (!is.null(sheets$weights)) study$weights <- matched_weights(sheets$weights, study$matrix)
    n_shared <- length(intersect(rownames(study$matrix), study$da$protein))
    ev_inform("{n_shared} of {nrow(study$matrix)} matrix proteins have DA results.",
      class = "enrichVolcano_study_overlap"
    )
  }
  structure(study, class = c("enrichVolcano_study", "list"))
}

read_sheets <- function(path) {
  if (dir.exists(path)) {
    files <- lapply(study_sheets, function(sheet) {
      found <- file.path(path, paste0(sheet, c(".csv", ".csv.gz")))
      found[file.exists(found)][1]
    })
    names(files) <- study_sheets
    files <- Filter(function(f) !is.na(f), files)
    return(lapply(files, utils::read.csv, check.names = FALSE))
  }
  if (!file.exists(path) || !grepl("\\.xlsx?$", path)) {
    ev_abort("{.path {path}} is not a study workbook or folder.", class = "enrichVolcano_input_error")
  }
  rlang::check_installed("readxl", reason = "to read study workbooks.")
  sheets <- intersect(study_sheets, readxl::excel_sheets(path))
  stats::setNames(lapply(sheets, function(s) as.data.frame(readxl::read_excel(path, s))), sheets)
}

check_samples <- function(samples) {
  require_columns(samples, c("sample", "group"))
  samples$sample <- as.character(samples$sample)
  twice <- unique(samples$sample[duplicated(samples$sample)])
  if (length(twice) > 0) {
    ev_abort("{.field samples} lists {.val {twice}} twice.", class = "enrichVolcano_input_error")
  }
  bad <- unique(samples$group[make.names(samples$group) != samples$group])
  if (length(bad) > 0) {
    ev_abort(
      "Group names must be syntactic to use in contrasts: {.val {bad}}.",
      class = "enrichVolcano_input_error"
    )
  }
  samples
}

check_contrasts <- function(contrasts, groups, da_contrasts) {
  require_columns(contrasts, c("name", "expression"))
  used <- unique(unlist(lapply(contrasts$expression, function(e) all.vars(str2lang(e)))))
  unknown <- setdiff(used, groups)
  if (length(unknown) > 0) {
    ev_abort("Contrasts use groups not in {.field samples}: {.val {unknown}}.", class = "enrichVolcano_input_error")
  }
  only_da <- setdiff(unique(da_contrasts), contrasts$name)
  only_sheet <- setdiff(contrasts$name, unique(da_contrasts))
  if (length(only_da) + length(only_sheet) > 0) {
    ev_abort(
      c(
        "{.field da_results} and {.field contrasts} must name the same contrasts.",
        i = "Only in da_results: {.val {only_da}}",
        i = "Only in contrasts: {.val {only_sheet}}"
      ),
      class = "enrichVolcano_input_error"
    )
  }
  contrasts
}

numeric_matrix <- function(tbl, name) {
  ids <- as.character(tbl[[1]])
  values <- tbl[-1]
  if (!all(vapply(values, is.numeric, logical(1)))) {
    ev_abort("Every sample column of {.field {name}} must be numeric.", class = "enrichVolcano_data_error")
  }
  if (anyDuplicated(ids) > 0) {
    ev_abort("{.field {name}} has duplicate protein IDs.", class = "enrichVolcano_data_error")
  }
  matrix(as.matrix(values), nrow(values), dimnames = list(ids, names(values)))
}

check_sample_links <- function(matrix, samples) {
  only_matrix <- setdiff(colnames(matrix), samples$sample)
  only_samples <- setdiff(samples$sample, colnames(matrix))
  if (length(only_matrix) + length(only_samples) > 0) {
    ev_abort(
      c(
        "Matrix columns and {.field samples$sample} must name the same samples.",
        i = "Only in matrix: {.val {only_matrix}}",
        i = "Only in samples: {.val {only_samples}}"
      ),
      class = "enrichVolcano_input_error"
    )
  }
}

matched_weights <- function(tbl, matrix) {
  weights <- numeric_matrix(tbl, "weights")
  same <- setequal(rownames(weights), rownames(matrix)) && setequal(colnames(weights), colnames(matrix)) &&
    identical(dim(weights), dim(matrix))
  if (!same) {
    ev_abort("{.field weights} must have the same proteins and samples as {.field matrix}.",
      class = "enrichVolcano_input_error"
    )
  }
  weights[rownames(matrix), colnames(matrix), drop = FALSE]
}

#' Example studies shipped with the package
#'
#' Seven proteomics studies from one lab, prepared as [read_study()] inputs.
#' The three limpa studies (`bfr_limpa`, `mouse_pas`, `yvo`) include the
#' matrix, samples,
#' contrasts and precision weights, so every test in [run_enrichment()] runs on
#' them; the rest ship DA results only, for fgsea and the figures. Call with no
#' name to list them with their species and designs.
#'
#' @param name A study name, or `NULL` to list them.
#' @return The study, as [read_study()] returns it, with its index row in
#'   `info`; or the index as a data frame.
#' @export
#' @examples
#' example_study()
example_study <- function(name = NULL) {
  root <- system.file("extdata", "studies", package = "enrichVolcano")
  index <- utils::read.csv(file.path(root, "index.csv"))
  if (is.null(name)) {
    return(index)
  }
  if (!isTRUE(name %in% index$name)) {
    ev_abort("No example study {.val {name}}; see {.code example_study()}.", class = "enrichVolcano_param_error")
  }
  info <- as.list(index[index$name == name, ])
  study <- read_study(file.path(root, name), species = info$species)
  study$info <- info
  study
}
