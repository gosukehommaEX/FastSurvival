#' Fast Kaplan-Meier Survival Probability at a Specified Time Point
#'
#' @description
#' Computes the Kaplan-Meier survival probability at a specified time point,
#' together with a standard error and confidence interval based on Greenwood's
#' variance formula. The function is a fully vectorized implementation that
#' avoids R-level loops and is designed for repeated use inside iterative
#' algorithms where many KM evaluations are required.
#'
#' @details
#' The KM estimate at time \code{t_eval} is defined as the product-limit
#' estimator evaluated at the largest observed event time less than or equal
#' to \code{t_eval}. If \code{t_eval} is smaller than the first observed
#' event time, \code{S(t) = 1} and the standard error is zero.
#'
#' The standard error is estimated by Greenwood's formula:
#'
#' SE[S(t)] = S(t) * sqrt(sum_{t_i <= t, d_i > 0} d_i / (n_i * (n_i - d_i)))
#'
#' where d_i is the number of events and n_i is the number at risk at time t_i.
#' The output field \code{std.err} follows the convention of
#' \code{\link[survival]{survfit}}, which reports SE[S(t)] / S(t) (i.e., the
#' standard error on the log scale) when \code{conf.type != "plain"}, and
#' SE[S(t)] when \code{conf.type = "plain"}. This function always returns
#' SE[S(t)] (the standard error on the survival scale).
#'
#' When \code{S(t_eval) = 0} (all subjects have experienced the event by
#' \code{t_eval}), the standard error is zero and the confidence interval
#' collapses to \code{[0, 0]}, consistent with \code{\link[survival]{survfit}}.
#'
#' When \code{presorted = TRUE} (default), the cutoff index is located via
#' \code{findInterval()}, which uses binary search and runs in O(log n).
#' When \code{presorted = FALSE}, the data are sorted internally before
#' computation; this is convenient for one-off calls but removes the
#' opportunity to reuse the same sorted arrays across many calls.
#'
#' Only observations at or before \code{t_eval} contribute to the KM product
#' and the Greenwood sum. Among those, only event positions (e_sorted == 1)
#' enter the Greenwood summation, so the inner loop over censored observations
#' is avoided entirely.
#'
#' Three confidence interval types are supported via \code{conf.type}:
#'
#' \itemize{
#'   \item \code{"plain"}: Linear interval on the survival scale,
#'     S(t) +/- z * SE. The bounds are clipped to [0, 1].
#'   \item \code{"log"}: Interval on the log scale (default in
#'     \code{\link[survival]{survfit}}),
#'     S(t) * exp(+/- z * SE / S(t)).
#'   \item \code{"log-log"}: Interval on the complementary log-log scale,
#'     S(t)^exp(+/- z * SE / (S(t) * log(S(t)))).
#' }
#'
#' @param t_sorted A numeric vector of event or censoring times. Must be
#'   sorted in ascending order when \code{presorted = TRUE}.
#' @param e_sorted An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{t_sorted}.
#' @param t_eval A single numeric value specifying the time point at which
#'   the survival probability is evaluated.
#' @param conf.int A single numeric value in (0, 1) specifying the
#'   confidence level. Defaults to 0.95.
#' @param conf.type A character string specifying the confidence interval type.
#'   Must be one of \code{"plain"}, \code{"log"}, or \code{"log-log"}.
#'   Defaults to \code{"log"}.
#' @param presorted A logical value. If \code{TRUE} (default), \code{t_sorted}
#'   and \code{e_sorted} are assumed to be sorted in ascending order of time,
#'   and \code{findInterval()} is used for O(log n) index lookup. If
#'   \code{FALSE}, the vectors are sorted internally before computation.
#'
#' @return A named numeric vector of length 4 with elements \code{surv},
#'   \code{std.err}, \code{lower}, and \code{upper}, representing the KM
#'   survival estimate, the Greenwood standard error SE[S(t)], and the lower
#'   and upper confidence limits at \code{t_eval}.
#'   Returns \code{c(surv = NA_real_, std.err = NA_real_, lower = NA_real_,
#'   upper = NA_real_)} when \code{n} is zero or \code{t_eval} precedes all
#'   observations.
#'
#' @examples
#' set.seed(42)
#' t_raw <- rexp(100, rate = 1 / 10)
#' e_raw <- rbinom(100, 1, 0.7)
#'
#' # presorted = TRUE (default): sort once outside, reuse inside a loop
#' ord <- order(t_raw)
#' t_s <- t_raw[ord]
#' e_s <- e_raw[ord]
#' survfit_fast(t_s, e_s, t_eval = 10, conf.type = "plain")
#' survfit_fast(t_s, e_s, t_eval = 10, conf.type = "log")
#' survfit_fast(t_s, e_s, t_eval = 10, conf.type = "log-log")
#'
#' # presorted = FALSE: sort internally, convenient for one-off calls
#' survfit_fast(t_raw, e_raw, t_eval = 10, presorted = FALSE)
#'
#' # Validation against survival::survfit
#' library(survival)
#' fit <- survfit(Surv(t_raw, e_raw) ~ 1, conf.type = "plain")
#' summary(fit, times = 10)
#'
#' @seealso
#' \code{\link[survival]{survfit}} for the standard KM estimator.
#'
#' @export
survfit_fast <- function(t_sorted, e_sorted, t_eval,
                         conf.int = 0.95, conf.type = "log",
                         presorted = TRUE) {
  conf.type <- match.arg(conf.type, choices = c("plain", "log", "log-log"))
  na_out <- c(surv = NA_real_, std.err = NA_real_,
              lower = NA_real_, upper = NA_real_)

  n <- length(t_sorted)
  if (n == 0L) return(na_out)

  # Sort internally when presorted = FALSE
  if (!presorted) {
    ord      <- order(t_sorted)
    t_sorted <- t_sorted[ord]
    e_sorted <- e_sorted[ord]
  }

  # Binary search for the last index <= t_eval (O(log n) when presorted)
  m <- findInterval(t_eval, t_sorted)
  if (m == 0L) {
    # t_eval is before any observation: S(t) = 1 by convention
    return(c(surv = 1, std.err = 0, lower = 1, upper = 1))
  }

  # Restrict to observations up to t_eval
  e_m <- e_sorted[seq_len(m)]

  # Event positions within 1..m (avoids processing censored rows)
  ev_pos  <- which(e_m == 1L)
  n_ev_m  <- length(ev_pos)

  # At-risk count at each event position: n - (pos - 1)
  n_risk_ev <- n - ev_pos + 1L

  # KM survival estimate: product only over event positions
  surv <- prod(1 - 1L / n_risk_ev)

  # S(t) = 0: standard error is zero, CI collapses to [0, 0]
  if (surv == 0) return(c(surv = 0, std.err = 0, lower = 0, upper = 0))

  # Greenwood standard error: SE[S(t)] = S(t) * sqrt(sum 1 / (n_i * (n_i - 1)))
  gw_sum  <- sum(1 / (n_risk_ev * (n_risk_ev - 1L)))
  std.err <- surv * sqrt(gw_sum)

  z <- qnorm(1 - (1 - conf.int) / 2)

  # Confidence interval by type
  ci <- if (conf.type == "plain") {
    c(lower = max(0, surv - z * std.err),
      upper = min(1, surv + z * std.err))
  } else if (conf.type == "log") {
    c(lower = surv * exp(-z * std.err / surv),
      upper = surv * exp( z * std.err / surv))
  } else {
    # log-log scale
    if (surv >= 1) return(c(surv = surv, std.err = std.err,
                            lower = NA_real_, upper = NA_real_))
    log_s    <- log(surv)
    se_ll    <- std.err / (surv * abs(log_s))
    theta_lo <- log(-log_s) - z * se_ll
    theta_hi <- log(-log_s) + z * se_ll
    c(lower = exp(-exp(theta_hi)),
      upper = exp(-exp(theta_lo)))
  }

  c(surv = surv, std.err = std.err, ci)
}
