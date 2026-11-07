#' Fast Cox Proportional Hazards Estimation for Two Groups
#'
#' @description
#' Fits the Cox proportional hazards model with a single treatment indicator
#' to a two-group parallel trial and returns the maximum partial likelihood
#' estimate of the hazard ratio, the standard error of its logarithm, and a
#' Wald-type confidence interval, using output names consistent with
#' \code{summary(survival::coxph(...))}. The C++ backend forms the risk sets
#' in a single pass over the pooled sorted data and runs the Newton-Raphson
#' iteration on the numbers at risk and of events at the distinct event
#' times, so that an iteration takes time in proportion to the number of
#' distinct event times rather than the number of subjects.
#'
#' @details
#' Let t_k (k = 1, ..., K) denote the distinct observed event times in the
#' pooled sample. At each t_k, let n_T_k and n_C_k be the numbers at risk in
#' the treatment and control groups just before t_k (subjects censored at t_k
#' are at risk at t_k), and let O_T_k, O_C_k, and O_k = O_T_k + O_C_k be the
#' numbers of events in the treatment group, in the control group, and in
#' total at t_k. With the log hazard ratio beta and r = exp(beta), the log
#' partial likelihood of Cox (1972) with the Breslow (1974) handling of tied
#' event times is
#'
#' l(beta) = sum_k [O_T_k beta - O_k log(n_C_k + n_T_k r)].
#'
#' The Efron (1977) handling treats the O_k tied events at t_k as occurring
#' one after another and replaces the k-th term by O_T_k beta minus the sum,
#' over m = 0, ..., O_k - 1, of log(n_C_k + n_T_k r - (m / O_k) (O_C_k + O_T_k
#' r)). Without tied event times the two coincide. With \code{strata}, the
#' risk sets are formed within each stratum and l(beta) is the sum over
#' strata, so that the hazard ratio is common to all strata and each stratum
#' has its own baseline hazard.
#'
#' The estimate maximizes l(beta) by Newton-Raphson iteration from beta = 0,
#' which Therneau and Grambsch (2000, Appendix E) work through on a small
#' data set. A step after which l(beta) decreases is halved. From the first
#' halved step on, the iteration also keeps the interval between the largest
#' value at which the score l'(beta) has been found positive and the smallest
#' value at which it has been found negative, and a trial value outside this
#' interval is replaced by its midpoint; this keeps the iteration convergent
#' when the hazard ratio is extreme and one group is small. The iteration
#' stops when the relative change in l(beta) is at most 1e-9, and at most 20
#' iterations are taken. The starting value, the convergence criterion, and
#' the maximum number of iterations are those of \code{survival::coxph()}
#' with the default \code{\link[survival]{coxph.control}()}, so that the two
#' estimates agree to rounding error when no step is halved; when a step is
#' halved, they can differ within the convergence tolerance.
#'
#' The standard error of the log hazard ratio is 1 / sqrt(I), where
#' I = -l''(beta) is the observed information at the estimate, and the Wald
#' confidence interval for the hazard ratio is exp(coef -/+ z se), where z is
#' the standard normal quantile for \code{conf.level}.
#'
#' Let O_T be the total number of events in the treatment group, L the number
#' of events at times when no control subject is at risk, and H the number of
#' events at times when a treatment-group subject is at risk. The score
#' decreases in beta, from O_T - L as beta tends to -Inf to O_T - H as beta
#' tends to Inf, so the estimate is finite if and only if L < O_T < H.
#' Otherwise (for example, when all events are in one group) the estimate is
#' infinite and the result is \code{NA}, whereas \code{survival::coxph()}
#' returns a large finite value, usually with a warning.
#'
#' Times are tied only when they are exactly equal. By default
#' \code{survival::coxph()} also treats as tied the times that differ only by
#' rounding error (\code{timefix = TRUE} in
#' \code{\link[survival]{coxph.control}()}); the two give the same result
#' when no distinct times are that close, and
#' \code{coxph.control(timefix = FALSE)} gives the comparable fit otherwise.
#'
#' The returned object has class \code{"coxph_fast"} and is a named numeric
#' vector of length 5. A \code{print()} method formats the result similarly
#' to \code{summary(coxph(...))}.
#'
#' @param time A numeric vector of follow-up times for all subjects (pooled
#'   over both groups).
#'   Alternatively, a formula \code{Surv(time, event) ~ group}, optionally with
#'   \code{strata()} terms, with the variables in
#'   \code{data}; see \code{\link{formula_interface}}.
#' @param event An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{time}.
#' @param group A vector of group labels aligned with \code{time}. Any type
#'   that supports equality comparison is accepted.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group. Subjects with \code{group != control} are
#'   treated as the treatment group.
#'   With the formula interface it defaults to the first level of the
#'   group variable.
#' @param side 1 for a one-sided test in the direction of treatment benefit
#'   (hazard ratio below 1, i.e. a negative coefficient) or 2 for a two-sided
#'   test (default 2). The returned vector has no p-value: \code{side} is
#'   stored as an attribute, and the \code{print()} method reports the
#'   p-value that follows it. The confidence interval is always two-sided at
#'   \code{conf.level}.
#' @param conf.level A single numeric value in (0, 1) specifying the confidence
#'   level for the Wald interval. Defaults to 0.95.
#' @param presorted A logical value. If \code{TRUE}, \code{time},
#'   \code{event}, and \code{group} are assumed to be already sorted in
#'   ascending order of \code{time} (with \code{strata}, sorted by stratum
#'   and by time within stratum, so that the rows of each stratum are
#'   contiguous), and the internal \code{order()} call is skipped. If
#'   \code{FALSE} (default), sorting is handled internally.
#'   The order is checked, and an error is given when it does not hold.
#' @param strata An optional vector of stratum labels aligned with
#'   \code{time}, without missing values. When supplied, the stratified
#'   estimator described in Details is computed. Several stratification
#'   factors can be combined with \code{interaction()}. \code{NULL} (default)
#'   gives the unstratified estimator.
#' @param data A data frame in which the variables of a formula given as
#'   \code{time} are evaluated. Used only with the formula interface; see
#'   \code{\link{formula_interface}}.
#' @param ties A character string naming the handling of tied event times,
#'   \code{"efron"} (default, as in \code{survival::coxph()}) or
#'   \code{"breslow"}. Without tied event times the two give the same result.
#'
#' @return An object of class \code{"coxph_fast"}, which is a named numeric
#'   vector of length 5 with elements matching the column names of
#'   \code{summary(coxph(...))$coefficients} and
#'   \code{summary(coxph(...))$conf.int}:
#' \describe{
#'   \item{\code{coef}}{Estimate of the log hazard ratio.}
#'   \item{\code{exp(coef)}}{Estimate of the hazard ratio.}
#'   \item{\code{se(coef)}}{Standard error of \code{coef}, from the observed
#'     information at the estimate.}
#'   \item{\code{lower .95}}{Lower bound of the Wald confidence interval for
#'     the hazard ratio. The label reflects \code{conf.level} (e.g.,
#'     \code{"lower .90"} when \code{conf.level = 0.90}).}
#'   \item{\code{upper .95}}{Upper bound of the Wald confidence interval.}
#' }
#' The attributes \code{ties}, \code{iter} (the number of Newton-Raphson
#' iterations), and \code{converged} record the fit. With \code{strata}, the
#' number of strata is stored in the attribute \code{strata}.
#' The elements are \code{NA_real_} (still with class \code{"coxph_fast"})
#' when the estimate cannot be computed (no events, or an infinite estimate
#' as described in Details), with \code{converged = NA}, and when the
#' iteration does not converge in 20 iterations, with \code{converged =
#' FALSE} and a warning.
#'   The object also has the class \code{"fastsurv"}, whose
#'   \code{as.data.frame()}, \code{coef()}, and \code{confint()} methods
#'   are described in \code{\link{as.data.frame.fastsurv}}.
#'
#' @examples
#' library(survival)
#'
#' # Compare coxph_fast with coxph on the ovarian dataset.
#' # coxph() treats rx as numeric with rx = 1 as the reference (control),
#' # so set control = 1 for a consistent comparison.
#' fit_fast <- coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
#' fit_fast
#'
#' fit_cox <- coxph(Surv(futime, fustat) ~ rx, data = ovarian)
#' c(coxph_fast = unname(fit_fast["coef"]), coxph = unname(coef(fit_cox)))
#' c(coxph_fast = unname(fit_fast["se(coef)"]),
#'   coxph = sqrt(fit_cox$var[1, 1]))
#'
#' # The Breslow handling of tied event times
#' coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1,
#'            ties = "breslow")
#'
#' # Stratified by residual disease, compared with coxph(... + strata())
#' coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1,
#'            strata = ovarian$resid.ds)
#' coef(coxph(Surv(futime, fustat) ~ rx + strata(resid.ds), data = ovarian))
#'
#' # presorted = TRUE: sort once outside, reuse inside a loop
#' ord <- order(ovarian$futime)
#' coxph_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
#'            control = 1, presorted = TRUE)
#'
#' # Formula interface, and the result as a data frame (log hazard ratio)
#' fit_f <- coxph_fast(Surv(futime, fustat) ~ rx, data = ovarian)
#' as.data.frame(fit_f)
#' exp(confint(fit_f))
#'
#' \donttest{
#' # Speed comparison against coxph()
#' if (requireNamespace("microbenchmark", quietly = TRUE)) {
#'   microbenchmark::microbenchmark(
#'     coxph_fast = coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 1),
#'     coxph      = coxph(Surv(futime, fustat) ~ rx, data = ovarian),
#'     times = 1000
#'   )
#' }
#' }
#'
#' @seealso
#' \code{\link[survival]{coxph}} for the Cox model with any covariates.
#' \code{\link{print.coxph_fast}} for the print method.
#'
#' @references
#' Breslow, N. (1974). Covariance analysis of censored survival data.
#' \emph{Biometrics}, \emph{30}(1), 89-99.
#'
#' Cox, D. R. (1972). Regression models and life-tables. \emph{Journal of the
#' Royal Statistical Society. Series B (Methodological)}, \emph{34}(2),
#' 187-202.
#'
#' Efron, B. (1977). The efficiency of Cox's likelihood function for censored
#' data. \emph{Journal of the American Statistical Association}, \emph{72}(359),
#' 557-565.
#'
#' Therneau, T. M., & Grambsch, P. M. (2000). \emph{Modeling Survival Data:
#' Extending the Cox Model}. New York: Springer.
#'
#' @importFrom stats qnorm setNames
#' @export
coxph_fast <- function(time, event, group, control, side = 2,
                       conf.level = 0.95, presorted = FALSE, strata = NULL,
                       data = NULL, ties = c("efron", "breslow")) {

  # Formula interface: Surv(time, event) ~ group + strata(s) evaluated in
  # 'data'
  if (inherits(time, "formula")) {
    f <- fastsurv_formula(time, data,
                          event = if (!missing(event)) event,
                          group = if (!missing(group)) group,
                          strata = strata, group_mode = "required",
                          strata_ok = TRUE, fun = "coxph_fast()")
    if (missing(control) || is.null(control)) control <- f$control
    return(coxph_fast(f$time, f$event, f$group, control = control,
                      side = side, conf.level = conf.level,
                      presorted = presorted, strata = f$strata,
                      ties = ties))
  }
  if (!is.null(data)) stop("'data' is used only with a formula")
  ties <- match.arg(ties)

  if (!side %in% c(1L, 2L)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided)")
  }
  if (length(conf.level) != 1L || !is.finite(conf.level) ||
      conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be in (0, 1)")
  }

  # Prepare NA output with coxph-compatible names
  ci_lab <- conf.level * 100
  ci_lo  <- sprintf("lower .%g", ci_lab)
  ci_hi  <- sprintf("upper .%g", ci_lab)
  na_out <- setNames(
    rep(NA_real_, 5L),
    c("coef", "exp(coef)", "se(coef)", ci_lo, ci_hi)
  )

  # Input validation
  n <- length(time)
  if (length(event) != n || length(group) != n) {
    stop("'time', 'event', and 'group' must have the same length")
  }
  use_strata <- !is.null(strata)
  if (use_strata) {
    if (length(strata) != n) {
      stop("'strata' must have the same length as 'time'")
    }
    if (anyNA(strata)) {
      stop("'strata' must not contain missing values")
    }
  }

  # Attach the attributes of a result (the number of strata when stratified)
  wrap <- function(v, n_strata = NULL, iter = 0L, converged = NA) {
    out <- structure(v, conf.level = conf.level, side = side,
                     control = control, ties = ties, iter = iter,
                     converged = converged,
                     class = c("coxph_fast", "fastsurv"))
    if (use_strata) attr(out, "strata") <- n_strata
    out
  }

  if (n == 0L) return(wrap(na_out, 0L))
  check_time_event(time, event)

  # Treatment indicator: 1 = treatment, 0 = control
  j <- two_group_indicator(group, control)

  # Stratum labels mapped to contiguous integers 1..S
  if (use_strata) {
    strata_int <- match(strata, sort(unique(strata)))
    n_strata   <- max(strata_int)
  } else {
    n_strata <- NULL
  }

  if (sum(event) == 0L) return(wrap(na_out, n_strata))

  # Sort by time (by stratum, then time, when stratified) when not presorted
  if (!presorted) {
    ord   <- if (use_strata) order(strata_int, time) else order(time)
    time  <- time[ord]
    event <- as.integer(event[ord])
    j     <- j[ord]
    if (use_strata) strata_int <- strata_int[ord]
  } else {
    check_presorted(time, if (use_strata) strata_int else NULL)
    event <- as.integer(event)
  }

  # C++ core: one scan over the pooled sorted data (per stratum when
  # stratified) and Newton-Raphson iteration on the counts at the event times
  # -> c(coef, se, loglik at 0, loglik at coef, iterations, status)
  res <- cox_core(time, event, j,
                  if (use_strata) strata_int else integer(0),
                  efron = ties == "efron")
  iter <- as.integer(res[5L])

  # Status 2: infinite or undefined estimate; status 1: not converged
  if (res[6L] == 2) return(wrap(na_out, n_strata, iter, NA))
  if (res[6L] == 1) {
    warning("the Newton-Raphson iteration did not converge in ", iter,
            " iterations; NA is returned", call. = FALSE)
    return(wrap(na_out, n_strata, iter, FALSE))
  }

  # Wald SE and CI on the log scale
  coef    <- res[1L]
  se_coef <- res[2L]
  z       <- qnorm(1 - (1 - conf.level) / 2)

  out <- setNames(
    c(coef, exp(coef), se_coef,
      exp(coef - z * se_coef),
      exp(coef + z * se_coef)),
    c("coef", "exp(coef)", "se(coef)", ci_lo, ci_hi)
  )

  wrap(out, n_strata, iter, TRUE)
}
