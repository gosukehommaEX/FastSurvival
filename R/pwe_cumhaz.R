#' Cumulative hazard of a piecewise-constant hazard
#'
#' Cumulative hazard at 'times' for the same parametrization. Built so that the
#' break-point cumulative values equal simdata_fast's piecewise_precompute().
#'
#' @param times Numeric vector of times.
#' @param haz Hazard vector of length K; a single value gives a constant hazard.
#' @param etime Breakpoints of length K + 1, from 0 to \code{Inf}; not used for
#'   a single hazard.
#' @return Numeric vector of the cumulative hazards at \code{times}.
#' @noRd
pwe_cumhaz <- function(times, haz, etime) {
  if (length(haz) == 1L) return(as.numeric(haz) * pmax(times, 0))
  k_int       <- length(haz)
  brk         <- etime[seq_len(k_int)]
  seglen      <- diff(brk)
  cumh_at_brk <- c(0, cumsum(haz[seq_len(k_int - 1L)] * seglen))
  idx         <- findInterval(times, brk)
  idx[idx < 1L] <- 1L
  cumh_at_brk[idx] + haz[idx] * (times - brk[idx])
}
