#' Launch the interactive enrichVolcano GUI
#'
#' Opens a Shiny app to load differential abundance results, choose databases
#' and thresholds, and render the volcano-in-ring composite interactively.
#' Requires the `shiny` package (Suggests).
#'
#' @param ... Passed to [shiny::runApp()] (e.g. `launch.browser`, `port`).
#' @return Invisibly `NULL`; called for the side effect of running the app.
#' @export
ev_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    ev_abort(c("The {.pkg shiny} package is required for {.fn ev_app}.",
               "i" = "Install it with {.code install.packages(\"shiny\")}."),
             class = "ev_shiny_missing")
  }
  app_dir <- system.file("shiny", package = "enrichVolcano")
  if (!nzchar(app_dir)) {
    ev_abort("Could not locate the bundled Shiny app.", class = "ev_app_missing")
  }
  shiny::runApp(app_dir, ...)
}
