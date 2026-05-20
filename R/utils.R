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
