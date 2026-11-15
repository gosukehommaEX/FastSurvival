#' Resolve whole-trial accrual for the multi-arm path
#'
#' Completes the accrual breakpoints (solving an open final interval from the
#' total when a trailing 'a.rate' is supplied) and returns the per-interval
#' accrual proportions, so every arm can be generated over a common accrual
#' window by passing the completed breakpoints together with these proportions
#' as a fully specified 'a.prop'. This mirrors the accrual resolution of the
#' single-endpoint path and is deliberately self-contained so the two-group path
#' is left byte-identical.
#'
#' @param n_total Total sample size.
#' @param a.time,a.rate,a.prop Accrual arguments of \code{simdata_fast()}.
#' @return A list with the completed breakpoints \code{a.time_full} and the
#'   proportions of the intervals \code{a_int_prob}.
#' @noRd
resolve_karm_accrual <- function(n_total, a.time, a.rate, a.prop) {
  n_int_time <- length(a.time) - 1L
  acc_tol    <- 1e-8 * max(1, n_total)
  use_rate   <- !is.null(a.rate)
  use_prop   <- !is.null(a.prop)
  if (use_rate == use_prop) stop("Supply exactly one of 'a.rate' and 'a.prop'")
  if (length(a.time) < 2L) stop("'a.time' must have at least two elements")
  if (any(diff(a.time) <= 0)) stop("'a.time' must be strictly increasing")

  if (use_rate) {
    if (any(a.rate <= 0)) stop("All 'a.rate' values must be positive")
    if (length(a.rate) == n_int_time) {
      implied <- sum(a.rate * diff(a.time))
      if (abs(implied - n_total) > acc_tol) {
        stop("'a.rate' implies ", round(implied, 4), " subjects over the ",
             "accrual period but sum(n) is ", n_total, ".")
      }
      a.time_full <- a.time
    } else if (length(a.rate) == n_int_time + 1L) {
      bounded <- if (n_int_time >= 1L) {
        sum(a.rate[seq_len(n_int_time)] * diff(a.time))
      } else {
        0
      }
      remaining <- n_total - bounded
      if (remaining <= 0) {
        stop("The specified accrual intervals already accrue ", round(bounded, 4),
             " subjects, at least sum(n) = ", n_total, ".")
      }
      final_dur   <- remaining / a.rate[length(a.rate)]
      a.time_full <- c(a.time, a.time[length(a.time)] + final_dur)
    } else {
      stop("'a.rate' must have length equal to length(a.time) - 1 or length(a.time)")
    }
    weights_a <- a.rate * diff(a.time_full)
  } else {
    if (any(a.prop <= 0)) stop("All 'a.prop' values must be positive")
    if (length(a.prop) != n_int_time) {
      stop("'a.prop' must have length equal to length(a.time) - 1")
    }
    a.time_full <- a.time
    weights_a   <- as.numeric(a.prop)
  }
  list(a.time_full = a.time_full, a_int_prob = weights_a / sum(weights_a))
}
