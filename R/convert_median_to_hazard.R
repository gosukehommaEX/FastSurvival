#' Convert median survival time to hazard
#'
#' @param median_arg Medians: a vector, or a (nested) list of vectors.
#' @return The hazards \code{log(2) / median}, with the same structure.
#' @noRd
convert_median_to_hazard <- function(median_arg) {
  if (is.list(median_arg)) {
    lapply(median_arg, convert_median_to_hazard)
  } else {
    log(2) / median_arg
  }
}
