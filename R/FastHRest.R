#' Calculate a hazard ratio for a two-arm parallel group trial.
#'
#' This function calculates the hazard ratio by various methods.
#'
#' @param time A numeric vector representing the event time.
#' @param event A numeric vector representing the flag of event for both groups.
#' @param group A numeric vector representing the group indicator.
#' @param control A numeric value what value/character represents the control group.
#' @param method Indicates the method: must be either "PY", "MST", "Pike", "Peto" or "LR".
#'
#' @return A dataframe includng a method and a hazard ratio.
#'
#' @examples
#' library(survival)
#' FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 1, 'PY')
#' FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 1, 'MST')
#' FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 1, 'Pike')
#' FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 1, 'Peto')
#' FastHRest(ovarian$futime, ovarian$fustat, ovarian$rx, 1, 'LR')
#'
#' @export
FastHRest = function(time, event, group, control, method) {
  # Convert groups to numeric values
  j = as.numeric(group != control)
  # Convert time to integer values
  time = match(time, sort(unique(time)))
  if(method == 'PY') {
    # PersonyYear method
    HR = (sum(event * j) / sum(time * j)) / ((sum(event * (1 - j)) / sum(time * (1 - j))))
  } else if(method == 'MST') {
    # Median survival method
    HR = median(time[j == 0]) / median(time[j == 1])
  } else {
    # Times events occurred
    t.k = sort(unique(time[event == 1]))
    # Define at risk
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
    # The observed number of events for the group j(=T,C) at t.k
    O1 = sum(e.1k)
    O0 = sum(e.0k)
    # The expected number of events for the group j at t.k
    E1 = sum(e.jk * (n.1k / n.jk))
    E0 = sum(e.jk * (n.0k / n.jk))
    # The variance
    V1 = sum((n.1k * n.0k * e.jk * (n.jk - e.jk)) / (n.jk ^ 2 * (n.jk - 1)), na.rm = TRUE)
    if(method == 'Pike') {
      # Pike method
      HR = (O1 * E0) / (O0 * E1)
    } else if(method == 'Peto') {
      # Peto method
      HR = exp((O1 - E1) / V1)
    } else if(method == 'LR') {
      # Log-rank test statistic
      LR = (O1 - E1) / sqrt(V1)
      # Log-rank test based method
      HR = exp(2 * LR / sqrt(O1 + O0))
    }
  }
  # Result
  result = data.frame(method, HR)
  return(result)
}
