#' Resolve accrual interval counts (illness-death path)
#'
#' Replicates the accrual resolution of the single-endpoint path for the
#' illness-death wrapper, returning the completed accrual breakpoints and the
#' per-interval subject counts for each group.
#'
#' @param n_grp_int Sample sizes of the groups (length 1 or 2).
#' @param a.time,a.rate,a.prop Accrual arguments of \code{simdata_fast()}.
#' @return A list with the completed breakpoints \code{a.time_full} and the
#'   counts per interval of the control and treatment groups,
#'   \code{acc_counts_c} and \code{acc_counts_t} (empty for one group).
#' @noRd
resolve_accrual_counts <- function(n_grp_int, a.time, a.rate, a.prop) {
  n_total    <- sum(n_grp_int)
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

  a_int_prob   <- weights_a / sum(weights_a)
  acc_counts_c <- accrual_cell_counts(as.integer(n_grp_int[1L]), a_int_prob)
  acc_counts_t <- if (length(n_grp_int) == 2L) {
    accrual_cell_counts(as.integer(n_grp_int[2L]), a_int_prob)
  } else {
    integer(0)
  }
  list(a.time_full = a.time_full,
       acc_counts_c = acc_counts_c, acc_counts_t = acc_counts_t)
}
