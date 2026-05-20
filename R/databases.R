#' List supported gene-set databases
#'
#' Returns a tibble with metadata for every registered database.
#'
#' @return Tibble with columns: name, display_name, type, species, source,
#'   license, citation_key, description_short.
#' @export
list_databases <- function() {
  path <- system.file("extdata", "db_registry.csv",
                      package = "enrichVolcano")
  if (!nzchar(path)) {
    ev_abort("Database registry not installed.", class = "ev_no_registry")
  }
  utils::read.csv(path, stringsAsFactors = FALSE) |>
    tibble::as_tibble()
}

#' Print extended description for one database
#'
#' @param name Database name (see `list_databases()$name`).
#' @return Invisibly, the registry row for `name`.
#' @export
database_info <- function(name) {
  reg <- list_databases()
  row <- reg[reg$name == name, ]
  if (nrow(row) == 0) {
    ev_abort("Unknown database {.val {name}}. See {.code list_databases()}.",
             class = "ev_unknown_database")
  }
  cli::cli_h2("{row$display_name} ({row$name})")
  cli::cli_dl(c(
    Type        = row$type,
    Species     = row$species,
    Source      = row$source,
    License     = row$license,
    Description = row$description_short
  ))
  invisible(row)
}
