#' @useDynLib sumgamma, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("sumgamma", libpath)
}
