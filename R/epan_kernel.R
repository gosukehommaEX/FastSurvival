#' Epanechnikov kernel
#'
#' @param u Numeric vector.
#' @return The kernel 0.75 (1 - u^2) for |u| <= 1, and 0 otherwise.
#' @noRd
epan_kernel <- function(u) ifelse(abs(u) <= 1, 0.75 * (1 - u * u), 0)
