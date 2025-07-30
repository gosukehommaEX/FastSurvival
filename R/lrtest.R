#' Calculate log-rank test statistic for a two-arm parallel group trial
#'
#' This function calculates the log-rank test statistic for comparing survival curves
#' between two groups. It can return either one-sided or two-sided test statistics.
#' The implementation is optimized for performance compared to standard survival
#' analysis functions, with optional tie handling for maximum speed.
#'
#' @param time A numeric vector representing the event times for all subjects.
#' @param event A numeric vector representing the event indicator for all subjects
#'   (1 = event occurred, 0 = censored).
#' @param group A numeric vector representing the group indicator for all subjects.
#' @param control A numeric value indicating which value in \code{group} represents
#'   the control group.
#' @param side A numeric value indicating the type of test: 1 for one-sided test,
#'   2 for two-sided test (returns chi-square statistic).
#' @param handling_tie A logical value indicating whether to handle tied observations.
#'   Set to FALSE for simulation data with continuous distributions (faster). Default is TRUE.
#' @param is_sorted A logical value indicating whether the data is already sorted by time.
#'   Set to TRUE for pre-sorted data to skip sorting operations (faster). Default is FALSE.
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
#' When \code{handling_tie = FALSE}, the function assumes no tied event times
#' (suitable for simulation data from continuous distributions) and uses a faster
#' algorithm that avoids ranking and tabulation operations.
#'
#' When \code{is_sorted = TRUE}, the function assumes the data is already sorted
#' by time and skips sorting operations for additional speed improvement.
#'
#' The function computes the test statistic using the standard log-rank formula:
#' \deqn{Z = \frac{O_1 - E_1}{\sqrt{V_1}}}
#' where \eqn{O_1} is the observed number of events in group 1, \eqn{E_1} is the
#' expected number of events in group 1, and \eqn{V_1} is the variance.
#'
#' @examples
#' library(survival)
#'
#' # One-sided log-rank test with tie handling
#' z_stat <- lrtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 1)
#' cat("One-sided Z statistic:", z_stat, "\n")
#'
#' # Two-sided log-rank test without tie handling (faster for simulation data)
#' set.seed(123)
#' n <- 1000
#' group_sim <- rep(c(1, 2), each = n/2)
#' time_sim <- c(rexp(n/2, 0.1), rexp(n/2, 0.15))
#' event_sim <- rep(1, n)
#' chisq_stat <- lrtest(time_sim, event_sim, group_sim, 1, 2, handling_tie = FALSE)
#' cat("Chi-square statistic:", chisq_stat, "\n")
#'
#' # Ultra-fast version for pre-sorted simulation data
#' ord <- order(time_sim)
#' chisq_stat_fast <- lrtest(time_sim[ord], event_sim[ord], group_sim[ord], 1, 2,
#'                          handling_tie = FALSE, is_sorted = TRUE)
#' cat("Ultra-fast result:", chisq_stat_fast, "\n")
#'
#' @seealso
#' \code{\link[survival]{survdiff}} for the standard implementation,
#' \code{\link{esthr}} for hazard ratio estimation
#'
#' @export
lrtest <- function(time, event, group, control, side, handling_tie = TRUE, is_sorted = FALSE) {
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

  # Convert groups to binary indicator (0 = control, 1 = treatment)
  j <- as.numeric(group != control)

  if (!handling_tie) {
    # Fast path for simulation data without ties
    if (is_sorted) {
      # Data is already sorted - use as is
      time_sorted <- time
      event_sorted <- event
      j_sorted <- j
    } else {
      # Sort by time
      ord <- order(time)
      time_sorted <- time[ord]
      event_sorted <- event[ord]
      j_sorted <- j[ord]
    }

    # Calculate observed and expected events
    n <- length(time)
    O1 <- sum(event_sorted * j_sorted)
    E1 <- sum(event_sorted * rev(cumsum(rev(j_sorted))) / (n:1))
    V1 <- sum(((n:1) - rev(cumsum(rev(j_sorted)))) * rev(cumsum(rev(j_sorted))) *
                event_sorted * ((n:1) - event_sorted) / ((n:1)^2 * ((n:1) - 1)),
              na.rm = TRUE)

  } else {
    # Standard path with tie handling
    if (is_sorted) {
      # Data is already sorted - rank without additional sorting
      time_ranked <- rank(time, ties.method = "min")
      # Get unique event times (already in sorted order)
      t.k <- unique(time_ranked[event == 1])
    } else {
      # Convert event times to ranks to handle tied observations
      time_ranked <- rank(time, ties.method = "min")
      # Get unique event times
      t.k <- sort.int(unique.default(time_ranked[event == 1]))
    }

    # Calculate numbers at risk at each event time
    # Using vectorized operations for efficiency
    n.1k <- rev(cumsum(rev(tabulate(time_ranked * j))))[t.k]
    n.1k[is.na(n.1k)] <- 0
    n.0k <- rev(cumsum(rev(tabulate(time_ranked * (1 - j)))))[t.k]
    n.0k[is.na(n.0k)] <- 0
    n.jk <- n.1k + n.0k

    # Number of events at each time point t.k
    e.1k <- tabulate(time_ranked * j * event)[t.k]
    e.1k[is.na(e.1k)] <- 0
    e.0k <- tabulate(time_ranked * (1 - j) * event)[t.k]
    e.0k[is.na(e.0k)] <- 0
    e.jk <- e.1k + e.0k

    # Calculate observed and expected events
    O1 <- sum(e.1k)
    E1 <- sum(e.jk * (n.1k / n.jk))

    # Calculate variance using optimized formula
    V1 <- sum((n.1k * n.0k * e.jk * (n.jk - e.jk)) / (n.jk^2 * (n.jk - 1)), na.rm = TRUE)
  }

  # Compute log-rank test statistic
  LR <- (O1 - E1) / sqrt(V1)

  # Return based on side argument
  if (side == 1) {
    return(LR)  # One-sided test statistic (Z-score)
  } else {
    return(LR^2)  # Two-sided test statistic (chi-square)
  }
}
