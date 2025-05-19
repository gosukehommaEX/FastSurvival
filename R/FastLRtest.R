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
  # Convert groups to numeric values 
  j = as.numeric(group != control)
  # Convert time to integer values
  time = match(time, sort(unique(time)))
  # Times occuring events
  t.k = sort(unique(time[event == 1]))
  # Set at risk
  n.1k = rev(cumsum(rev(tabulate(time * j))))[t.k]
  n.1k[is.na(n.1k)] = 0
  n.0k = rev(cumsum(rev(tabulate(time * (1 - j)))))[t.k]
  n.0k[is.na(n.0k)] = 0
  n.jk = n.1k + n.0k
  # Number of events at t.k
  e.1k = tabulate(time * j * event)[t.k]
  e.1k[is.na(e.1k)] = 0
  e.0k = tabulate(time * (1 - j) * event)[t.k]
  e.0k[is.na(e.0k)] = 0
  e.jk = e.1k + e.0k
  # The observed number of events on the group 1 at t.k
  O1 = sum(e.1k)
  # The expected number of events on the group 1 at t.k
  E1 = sum(e.jk * (n.1k / n.jk))
  # The variance
  V1 = sum((n.1k * n.0k * e.jk * (n.jk - e.jk)) / (n.jk ^ 2 * (n.jk - 1)), na.rm = TRUE)
  # Log-ran.jk test statistic
  LR = (O1 - E1) / sqrt(V1)
  return(LR)
}
