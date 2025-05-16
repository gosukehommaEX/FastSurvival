#' Calculate a one-sided log-rank test statistic for a two-arm parallel group trial.
#'
#' This function calculates the one-sided log-rank test statistic.
#'
#' @param time A numeric vector representing the event time.
#' @param event A numeric vector representing the flag of event for both groups.
#' @param group A numeric vector representing the group indicator.
#' @param control A numeric value what value/character represents the control group.
#'
#' @return A numeric value of the one-sided log-rank test statistic.
#'
#' @examples
#' library(survival)
#' FastLRtest(ovarian$futime, ovarian$fustat, ovarian$rx, 2)
#'
#' @export
FastLRtest = function(time, event, group, control) {
  # Re-label treatment and control groups
  group = as.numeric(group != control)
  # Set the ordered distinct event times
  t.k = order(time)
  # Set the ordered events, groups, and sample size at t.k
  e.k = event[t.k]
  g.k = group[t.k]
  n.k = length(t.k):1
  n.1k = rev(cumsum(rev(g.k)))
  n.2k = n.k - n.1k
  # The observed number of events on the group 1 at t.k
  O1 = c(e.k %*% g.k)
  # The expected number of events on the group 1 at t.k
  E1 = c(e.k %*% (n.1k / n.k))
  # The variance
  V1 = sum(n.1k * n.2k * e.k * (n.k - e.k) / (n.k ^ 2 * (n.k - 1)), na.rm = TRUE)
  # Log-rank test statistic
  LR = (O1 - E1) / sqrt(V1)
  return(LR)
}
