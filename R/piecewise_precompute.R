#' Piecewise-exponential pre-computation
#'
#' Piecewise-exponential pre-computation: validates the hazard / breakpoint
#' pair and returns the hazard vector, the finite left breakpoints, and the
#' cumulative hazard at those breakpoints, matching rpiece_exp_r. For a single
#' hazard the breakpoint / cumulative entries are empty. 'nm' is the prefix of
#' the argument names used in the error messages ("e", "d", "h01", and so on).
#'
#' @param hazard Hazard vector.
#' @param e.time Breakpoints of length \code{length(hazard) + 1}, from 0 to
#'   \code{Inf}; not used for a single hazard.
#' @param nm Prefix of the argument names in the error messages.
#' @return A list with \code{hazard}, the finite left breakpoints
#'   \code{fin_time}, and the cumulative hazard \code{cum_haz} at them.
#' @noRd
piecewise_precompute <- function(hazard, e.time, nm = "e") {
  if (!is.numeric(hazard) || length(hazard) < 1L || anyNA(hazard) ||
      any(hazard < 0) || any(is.infinite(hazard))) {
    stop("Hazards must be finite and non-negative, and medians must be ",
         "positive")
  }
  if (length(hazard) == 1L) {
    return(list(hazard = as.numeric(hazard),
                fin_time = numeric(0), cum_haz = numeric(0)))
  }
  if (is.null(e.time)) {
    stop("'", nm, ".time' must be supplied when '", nm, ".hazard' (or '", nm,
         ".median') is a vector")
  }
  if (length(e.time) != length(hazard) + 1L) {
    stop("length(", nm, ".time) must equal length(", nm, ".hazard) + 1")
  }
  if (!is.infinite(e.time[length(e.time)])) {
    stop("Last element of '", nm, ".time' must be Inf")
  }
  if (anyNA(e.time) || e.time[1L] != 0 || any(diff(e.time) <= 0)) {
    stop("Piecewise breakpoints ('e.time', 'd.time', and so on) must start ",
         "at 0 and be strictly increasing")
  }
  n_int    <- length(hazard)
  fin_time <- e.time[seq_len(n_int)]
  lengths  <- diff(e.time[seq_len(n_int)])
  cum_haz  <- c(0, cumsum(hazard[-n_int] * lengths))
  list(hazard = as.numeric(hazard),
       fin_time = as.numeric(fin_time), cum_haz = as.numeric(cum_haz))
}
