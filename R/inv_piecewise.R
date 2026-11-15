#' Inverse cumulative hazard of a piecewise exponential
#'
#' 'pc' is the output of piecewise_precompute(); 'e' are unit exponential
#' draws. The time t solves H(t) = e, and is Inf when the last hazard is zero
#' and e exceeds the cumulative hazard at the last breakpoint.
#'
#' @param e Unit exponential draws.
#' @param pc The output of \code{piecewise_precompute()}.
#' @return The times t with H(t) = \code{e}; \code{Inf} where the cumulative
#'   hazard does not reach \code{e}.
#' @noRd
inv_piecewise <- function(e, pc) {
  if (length(pc$hazard) == 1L) {
    return(if (pc$hazard > 0) e / pc$hazard else rep(Inf, length(e)))
  }
  k <- findInterval(e, pc$cum_haz)
  h <- pc$hazard[k]
  ifelse(h > 0, pc$fin_time[k] + (e - pc$cum_haz[k]) / h, Inf)
}
