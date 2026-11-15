#' Median of a survival curve on a grid
#'
#' Median read off a survival curve by linear interpolation. Returns NA when the
#' curve does not reach 0.5 within the grid.
#'
#' @param times Time grid.
#' @param surv Survival values on the grid, or \code{NULL}.
#' @return The median, or \code{NA} when \code{surv} is \code{NULL} or does not
#'   reach 0.5.
#' @noRd
curve_median <- function(times, surv) {
  if (is.null(surv) || all(surv > 0.5)) return(NA_real_)
  i <- which(surv <= 0.5)[1L]
  if (i == 1L) return(times[1L])
  t0 <- times[i - 1L]; t1 <- times[i]
  s0 <- surv[i - 1L];  s1 <- surv[i]
  t0 + (0.5 - s0) * (t1 - t0) / (s1 - s0)
}
