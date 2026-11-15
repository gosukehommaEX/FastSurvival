#' Piecewise-constant hazard at given times
#'
#' Piecewise-constant hazard evaluated at 'times'. 'haz' is a length-K hazard
#' vector and 'etime' is the length-(K + 1) breakpoint vector with etime[1] == 0
#' and the last element Inf, matching simdata_fast. A single hazard is constant.
#'
#' @param times Numeric vector of times.
#' @param haz Hazard vector of length K; a single value gives a constant hazard.
#' @param etime Breakpoints of length K + 1, from 0 to \code{Inf}; not used for
#'   a single hazard.
#' @return Numeric vector of the hazards at \code{times}.
#' @noRd
pwe_hazard <- function(times, haz, etime) {
  if (length(haz) == 1L) return(rep(as.numeric(haz), length(times)))
  k_int <- length(haz)
  brk   <- etime[seq_len(k_int)]
  idx   <- findInterval(times, brk)
  idx[idx < 1L] <- 1L
  haz[idx]
}
