#' @keywords internal
#' @importFrom rlang .data %||%
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}
