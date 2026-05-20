.onLoad <- function(libname, pkgname) {
  options(
    enrichVolcano.cache_dir = tools::R_user_dir("enrichVolcano", which = "cache")
  )
  invisible()
}
