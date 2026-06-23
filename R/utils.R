#' Package-level imports
#'
#' @importFrom rlang .data
#' @importFrom rlang %||%
#' @keywords internal
#' @name enrichVolcano-imports
NULL

#' Abort with a classed cli error
#'
#' Forwards `.envir = rlang::caller_env()` so glue interpolation sees the
#' caller's locals (otherwise `{var}` resolves in `ev_abort`'s frame).
#' @keywords internal
#' @noRd
ev_abort <- function(message, ..., class = NULL,
                     call = rlang::caller_env(),
                     .envir = rlang::caller_env()) {
  cli::cli_abort(
    message,
    ...,
    class = c(class, "enrichVolcano_error"),
    call = call,
    .envir = .envir
  )
}

#' Inform with a classed cli message
#' @keywords internal
#' @noRd
ev_inform <- function(message, ..., class = NULL,
                      .envir = rlang::caller_env()) {
  cli::cli_inform(
    message,
    ...,
    class = c(class, "enrichVolcano_message"),
    .envir = .envir
  )
}

#' Warn with a classed cli warning
#' @keywords internal
#' @noRd
ev_warn <- function(message, ..., class = NULL,
                    .envir = rlang::caller_env()) {
  cli::cli_warn(
    message,
    ...,
    class = c(class, "enrichVolcano_warning"),
    .envir = .envir
  )
}

#' Standard abort for a missing column on a user data.frame
#'
#' @keywords internal
#' @noRd
ev_abort_missing_column <- function(df, col_name, col_arg, df_name) {
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  hints <- c(
    "i" = "Available columns: {.val {names(df)}}.",
    if (length(numeric_cols)) {
      c("i" = "Numeric columns: {.val {numeric_cols}}.")
    },
    "i" = "Pass {.code {col_arg} = \"<name>\"} or rename your column."
  )
  ev_abort(
    c(paste0("Column {.val {col_name}} not found in {.arg {df_name}}."), hints),
    class = "enrichVolcano_column_error"
  )
}

#' Assert that a value is a single recognisable colour
#'
#' `grDevices::col2rgb()` accepts names ("red"), hex ("#1B9E77"), and the
#' integer palette; anything else throws deep inside grid when the plot draws.
#' Catch it at the boundary instead.
#' @keywords internal
#' @noRd
ev_assert_colour <- function(x, arg = "disc_color") {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  ok <- is.character(x) && length(x) == 1L && !is.na(x) &&
    tryCatch(
      {
        grDevices::col2rgb(x)
        TRUE
      },
      error = function(e) FALSE
    )
  if (!ok) {
    ev_abort("{.arg {arg}} must be a single colour name or hex code; got {.val {x}}.",
      class = "ev_invalid_colour"
    )
  }
  invisible(x)
}
