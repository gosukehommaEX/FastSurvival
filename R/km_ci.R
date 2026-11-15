#' Pointwise confidence limits of a Kaplan-Meier curve
#'
#' Vectorized Kaplan-Meier confidence limits matching survfit_fast: "plain"
#' (survival scale), "log" (default), and "log-log" (complementary log-log).
#'
#' @param surv Kaplan-Meier estimates.
#' @param se Their standard errors on the survival scale.
#' @param conf.level Confidence level.
#' @param conf.type Type of the limits: \code{"plain"}, \code{"log"}, or
#'   \code{"log-log"}.
#' @return A list with the lower and upper limits, \code{lower} and
#'   \code{upper}, within 0 and 1.
#' @noRd
km_ci <- function(surv, se, conf.level, conf.type) {
  z     <- stats::qnorm(1 - (1 - conf.level) / 2)
  lower <- surv
  upper <- surv
  pos   <- surv > 0 & surv < 1 & is.finite(se)
  if (conf.type == "plain") {
    lower[pos] <- surv[pos] - z * se[pos]
    upper[pos] <- surv[pos] + z * se[pos]
  } else if (conf.type == "log") {
    lower[pos] <- surv[pos] * exp(-z * se[pos] / surv[pos])
    upper[pos] <- surv[pos] * exp( z * se[pos] / surv[pos])
  } else {
    s        <- surv[pos]
    log_s    <- log(s)
    se_ll    <- se[pos] / (s * abs(log_s))
    theta_lo <- log(-log_s) - z * se_ll
    theta_hi <- log(-log_s) + z * se_ll
    lower[pos] <- exp(-exp(theta_hi))
    upper[pos] <- exp(-exp(theta_lo))
  }
  zero <- surv <= 0
  lower[zero] <- 0
  upper[zero] <- 0
  lower <- pmax(0, pmin(1, lower))
  upper <- pmax(0, pmin(1, upper))
  list(lower = lower, upper = upper)
}
