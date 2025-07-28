#' Calculate log-rank test statistic for a two-arm parallel group trial
#'
#' This function calculates the log-rank test statistic for comparing survival curves
#' between two groups. It can return either one-sided or two-sided test statistics.
#'
#' @param time A numeric vector representing the event times for all subjects.
#' @param event A numeric vector representing the event indicator for all subjects
#'   (1 = event occurred, 0 = censored).
#' @param group A numeric vector representing the group indicator for all subjects.
#' @param control A numeric value indicating which value in \code{group} represents
#'   the control group.
#' @param side A numeric value indicating the type of test: 1 for one-sided test,
#'   2 for two-sided test (returns chi-square statistic).
#'
#' @return A numeric value of the log-rank test statistic. For \code{side = 1},
#'   returns the standardized log-rank statistic (Z-score). For \code{side = 2},
#'   returns the chi-square statistic (Z^2).
#'
#' @details
#' The log-rank test is a non-parametric test used to compare the survival
#' distributions of two groups. This implementation uses an optimized algorithm
#' that is significantly more efficient than the standard \code{\link[survival]{survdiff}}
#' function while producing equivalent results.
#'
#' The function converts event times to ranks to handle tied observations and
#' computes the test statistic using the standard log-rank formula:
#' \deqn{Z = \frac{O_1 - E_1}{\sqrt{V_1}}}
#' where \eqn{O_1} is the observed number of events in group 1, \eqn{E_1} is the
#' expected number of events in group 1, and \eqn{V_1} is the variance.
#'
#' @examples
#' library(survival)
#'
#' # One-sided log-rank test
#' z_stat <- lrtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 1)
#' cat("One-sided Z statistic:", z_stat, "\n")
#'
#' # Two-sided log-rank test (chi-square)
#' chisq_stat <- lrtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 2)
#' cat("Chi-square statistic:", chisq_stat, "\n")
#'
#' # Compare with survdiff
#' survdiff_result <- survdiff(Surv(futime, fustat) ~ rx, data = ovarian)
#' cat("survdiff chi-square:", survdiff_result$chisq, "\n")
#' cat("lrtest chi-square:", chisq_stat, "\n")
#'
#' # Performance comparison
#' \dontrun{
#' library(microbenchmark)
#' microbenchmark(
#'   survdiff = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
#'   lrtest = lrtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 2),
#'   times = 100
#' )
#' }
#'
#' @seealso
#' \code{\link[survival]{survdiff}} for the standard implementation,
#' \code{\link{esthr}} for hazard ratio estimation
#'
#' @export
lrtest <- function(time, event, group, control, side) {
  # Input validation
  if (length(time) != length(event) || length(time) != length(group)) {
    stop("Arguments 'time', 'event', and 'group' must have the same length")
  }

  if (!side %in% c(1, 2)) {
    stop("Argument 'side' must be either 1 (one-sided) or 2 (two-sided)")
  }

  if (sum(event) == 0) {
    stop("No events observed in the data")
  }

  # Convert groups to numeric values
  j <- as.numeric(group != control)

  # Convert time to integer values to handle ties
  time <- rank(time, ties.method = "min")

  # Times when events occurred
  t.k <- sort.int(unique.default(time[event == 1]))

  # Define numbers at risk
  n.1k <- rev(cumsum(rev(tabulate(time * j))))[t.k]
  n.1k[is.na(n.1k)] <- 0
  n.0k <- rev(cumsum(rev(tabulate(time * (1 - j)))))[t.k]
  n.0k[is.na(n.0k)] <- 0
  n.jk <- n.1k + n.0k

  # Number of events at each time point t.k
  e.1k <- tabulate(time * j * event)[t.k]
  e.1k[is.na(e.1k)] <- 0
  e.0k <- tabulate(time * (1 - j) * event)[t.k]
  e.0k[is.na(e.0k)] <- 0
  e.jk <- e.1k + e.0k

  # Calculate observed and expected events
  O1 <- sum(e.1k)
  O0 <- sum(e.0k)
  E1 <- sum(e.jk * (n.1k / n.jk))
  E0 <- sum(e.jk * (n.0k / n.jk))

  # Calculate variance
  V1 <- sum((n.1k * n.0k * e.jk * (n.jk - e.jk)) / (n.jk ^ 2 * (n.jk - 1)), na.rm = TRUE)

  # Log-rank test statistic
  LR <- (O1 - E1) / sqrt(V1)

  # Return based on side argument
  if (side == 1) {
    return(LR)  # One-sided test statistic (Z-score)
  } else if (side == 2) {
    return(LR^2)  # Two-sided test statistic (chi-square)
  }
}
