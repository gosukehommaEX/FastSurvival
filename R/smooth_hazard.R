#' Kernel-smoothed hazard from Nelson-Aalen increments
#'
#' @param grid Times at which the hazard is estimated.
#' @param te Event times.
#' @param dLambda Nelson-Aalen increments at \code{te}.
#' @param bw Bandwidth.
#' @return The hazard at \code{grid}; \code{NA} when there are no events.
#' @noRd
smooth_hazard <- function(grid, te, dLambda, bw) {
  if (length(te) == 0L) return(rep(NA_real_, length(grid)))
  vapply(grid, function(t0) {
    sum(epan_kernel((t0 - te) / bw) * dLambda) / bw
  }, numeric(1L))
}
