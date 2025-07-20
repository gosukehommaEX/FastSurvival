#' Calculate log-rank test statistic for a two-arm parallel group trial
#'
#' This function calculates the log-rank test statistic for comparing survival curves
#' between two groups. It can return either one-sided or two-sided test statistics.
#'
#' @param time A numeric vector representing the event time.
#' @param event A numeric vector representing the flag of event for both groups.
#' @param group A numeric vector representing the group indicator.
#' @param control A numeric value what value/character represents the control group.
#' @param side A numeric value indicating the type of test: 1 for one-sided, 2 for two-sided.
#'
#' @return A numeric value of the log-rank test statistic.
#'
#' @examples
#' library(survival)
#' # One-sided test
#' LRtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 1)
#' # Two-sided test
#' LRtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2, 2)
#'
#' @export
LRtest <- function(time, event, group, control, side) {
  # Convert groups to numeric values
  j <- as.numeric(group != control)

  # Convert time to integer values
  time <- rank(time, ties.method = "min")

  # Times events occurred
  t.k <- sort.int(unique.default(time[event == 1]))

  # Define at risk
  n.1k <- rev(cumsum(rev(tabulate(time * j))))[t.k]
  n.1k[is.na(n.1k)] <- 0
  n.0k <- rev(cumsum(rev(tabulate(time * (1 - j)))))[t.k]
  n.0k[is.na(n.0k)] <- 0
  n.jk <- n.1k + n.0k

  # Number of events at t.k
  e.1k <- tabulate(time * j * event)[t.k]
  e.1k[is.na(e.1k)] <- 0
  e.0k <- tabulate(time * (1 - j) * event)[t.k]
  e.0k[is.na(e.0k)] <- 0
  e.jk <- e.1k + e.0k

  # The observed number of events for the group j(=T,C) at t.k
  O1 <- sum(e.1k)
  O0 <- sum(e.0k)

  # The expected number of events for the group j at t.k
  E1 <- sum(e.jk * (n.1k / n.jk))
  E0 <- sum(e.jk * (n.0k / n.jk))

  # The variance
  V1 <- sum((n.1k * n.0k * e.jk * (n.jk - e.jk)) / (n.jk ^ 2 * (n.jk - 1)), na.rm = TRUE)

  # Log-rank test statistic
  LR <- (O1 - E1) / sqrt(V1)

  # Return based on side argument
  if (side == 1) {
    return(LR)  # One-sided test statistic
  } else if (side == 2) {
    return(LR^2)  # Two-sided test statistic (chi-square)
  }
}
