#' Compare Milestone Survival Probabilities Between Two Groups
#'
#' Compares the Kaplan-Meier survival probabilities of two groups at a
#' prespecified milestone timepoint. The point estimate of interest is the
#' difference in milestone survival, treatment minus control. Three inference
#' methods are provided. The \code{"wald"} method uses the unpooled Greenwood
#' variance directly. The \code{"loglog"} and \code{"mover"} methods build the
#' confidence interval for the difference with the method of variance estimates
#' recovery (MOVER), recovering the variance from the one-sample complementary
#' log-log and log transformed confidence intervals respectively. See Tang
#' (2021) for the MOVER difference interval and Tang (2022) for the use of
#' milestone survival in trial design.
#'
#' @param time A numeric vector of follow-up times.
#' @param status An integer vector of event indicators, 1 for an event and 0
#'   for a censored observation.
#' @param group A vector with exactly two distinct values identifying the
#'   group.
#' @param control The value of \code{group} that denotes the control group. The
#'   other value is the treatment group and the difference is reported as
#'   treatment minus control.
#' @param tau The milestone timepoint at which the survival probabilities are
#'   compared. A single positive number.
#' @param method The inference method for the difference in milestone survival,
#'   one of \code{"wald"}, \code{"loglog"}, or \code{"mover"}.
#' @param side The alternative hypothesis for the difference, one of
#'   \code{"two.sided"}, \code{"upper"} (treatment survival larger), or
#'   \code{"lower"} (treatment survival smaller). The confidence interval is
#'   always reported as a two-sided interval at \code{conf.level}.
#' @param conf.level The confidence level for the reported intervals.
#' @param presorted Logical. If \code{TRUE} the input is assumed to be sorted
#'   by \code{time} in ascending order and the internal sort is skipped. This
#'   is intended for repeated calls inside simulation loops.
#'
#' @return An object of class \code{"milestone_fast"}, a list with the
#'   per-group milestone survival estimates and standard errors, the difference
#'   estimate with its confidence interval, the test statistic, and the p-value.
#'
#' @references
#' Tang, Y. (2021). Some new confidence intervals for Kaplan-Meier based
#' estimators from one and two sample survival data. \emph{Statistics in
#' Medicine}, \emph{40}(23), 4961-4976.
#'
#' Tang, Y. (2022). Complex survival trial design by the product integration
#' method. \emph{Statistics in Medicine}, \emph{41}(4), 798-814.
#'
#' @examples
#' set.seed(1)
#' time <- c(rexp(50, 0.1), rexp(50, 0.07))
#' status <- rep(1, 100)
#' group <- rep(c(0, 1), each = 50)
#' milestone_fast(time, status, group, control = 0, tau = 10, method = "loglog")
#'
#' @export
milestone_fast <- function(time, status, group, control, tau,
                           method = c("wald", "loglog", "mover"),
                           side = c("two.sided", "upper", "lower"),
                           conf.level = 0.95, presorted = FALSE) {
  method <- match.arg(method)
  side <- match.arg(side)

  if (!is.numeric(time)) stop("'time' must be numeric.")
  if (length(status) != length(time) || length(group) != length(time)) {
    stop("'time', 'status', and 'group' must have the same length.")
  }
  status <- as.integer(status)
  if (any(is.na(status)) || any(!status %in% c(0L, 1L))) {
    stop("'status' must contain only 0 (censored) and 1 (event).")
  }
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("'tau' must be a single positive number.")
  }
  if (length(conf.level) != 1L || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number in (0, 1).")
  }

  # Map the two group levels to control (0) and treatment (1).
  if (is.factor(group)) {
    lev <- levels(droplevels(group))
  } else {
    lev <- sort(unique(group))
  }
  if (length(lev) != 2L) {
    stop("'group' must have exactly two distinct values.")
  }
  if (missing(control) || length(control) != 1L ||
      !(as.character(control) %in% as.character(lev))) {
    stop("'control' must be one of the two values in 'group'.")
  }
  trt <- lev[as.character(lev) != as.character(control)]
  grp01 <- as.integer(as.character(group) != as.character(control))

  core <- milestone_core(as.numeric(time), status, grp01,
                         as.numeric(tau), as.logical(presorted))

  surv0 <- core$surv0
  surv1 <- core$surv1
  v0 <- core$var0
  v1 <- core$var1
  se0 <- sqrt(v0)
  se1 <- sqrt(v1)

  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  # One-sample confidence limits for each group under the chosen transform.
  one_sample_ci <- function(surv, se) {
    if (is.na(surv) || is.na(se)) return(c(NA_real_, NA_real_))
    if (method == "wald") {
      return(c(surv - z * se, surv + z * se))
    }
    if (surv <= 0 || surv >= 1) return(c(NA_real_, NA_real_))
    if (method == "mover") {
      # log transform: var(log surv) = se^2 / surv^2
      g_se <- se / surv
      return(c(surv * exp(-z * g_se), surv * exp(z * g_se)))
    }
    # method == "loglog": g = log(-log(surv))
    g_se <- se / (surv * abs(log(surv)))
    return(c(surv^exp(z * g_se), surv^exp(-z * g_se)))
  }

  ci0 <- one_sample_ci(surv0, se0)
  ci1 <- one_sample_ci(surv1, se1)

  diff_est <- surv1 - surv0

  if (method == "wald") {
    se_diff <- sqrt(v0 + v1)
    diff_lower <- diff_est - z * se_diff
    diff_upper <- diff_est + z * se_diff
    stat <- diff_est / se_diff
  } else {
    # MOVER difference interval (Tang 2021, equation 10).
    l0 <- ci0[1]; u0 <- ci0[2]; l1 <- ci1[1]; u1 <- ci1[2]
    sigma_l <- (l1 - surv1)^2 + (u0 - surv0)^2
    sigma_u <- (u1 - surv1)^2 + (l0 - surv0)^2
    diff_lower <- diff_est - sqrt(sigma_l)
    diff_upper <- diff_est + sqrt(sigma_u)
    if (method == "loglog") {
      # Variance-stabilised complementary log-log test of equal survival.
      if (surv0 > 0 && surv0 < 1 && surv1 > 0 && surv1 < 1) {
        g0 <- log(-log(surv0))
        g1 <- log(-log(surv1))
        vg0 <- v0 / (surv0 * log(surv0))^2
        vg1 <- v1 / (surv1 * log(surv1))^2
        stat <- (g1 - g0) / sqrt(vg0 + vg1)
      } else {
        stat <- NA_real_
      }
    } else {
      # MOVER binding-side test statistic (confidence-level invariant).
      bind <- if (diff_est >= 0) sigma_l else sigma_u
      stat <- if (is.finite(bind) && bind > 0) z * diff_est / sqrt(bind) else NA_real_
    }
  }

  p_value <- if (is.na(stat)) {
    NA_real_
  } else if (side == "two.sided") {
    2 * stats::pnorm(-abs(stat))
  } else if (side == "upper") {
    stats::pnorm(-stat)
  } else {
    stats::pnorm(stat)
  }

  structure(
    list(
      surv = c(control = surv0, treatment = surv1),
      std.err = c(control = se0, treatment = se1),
      surv.lower = c(control = ci0[1], treatment = ci1[1]),
      surv.upper = c(control = ci0[2], treatment = ci1[2]),
      n = c(control = core$n0, treatment = core$n1),
      n.event = c(control = core$n_event0, treatment = core$n_event1),
      n.risk = c(control = core$n_risk0, treatment = core$n_risk1),
      diff = diff_est,
      diff.lower = diff_lower,
      diff.upper = diff_upper,
      statistic = stat,
      p.value = p_value,
      tau = as.numeric(tau),
      method = method,
      side = side,
      conf.level = conf.level,
      group.labels = c(control = as.character(control),
                       treatment = as.character(trt))
    ),
    class = "milestone_fast"
  )
}
