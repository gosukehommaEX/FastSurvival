#' Estimate hazard ratio for a two-arm parallel group trial
#'
#' This function estimates the hazard ratio using various established methods.
#' It provides multiple approaches for hazard ratio estimation, each with different
#' theoretical foundations and computational characteristics.
#'
#' @param time A numeric vector representing the event times for all subjects.
#' @param event A numeric vector representing the event indicator for all subjects
#'   (1 = event occurred, 0 = censored).
#' @param group A numeric vector representing the group indicator for all subjects.
#' @param control A numeric value indicating which value in \code{group} represents
#'   the control group.
#' @param method A character string specifying the estimation method. Must be one of:
#'   \itemize{
#'     \item \code{"PY"}: Person-Year method
#'     \item \code{"Pike"}: Pike's method
#'     \item \code{"Peto"}: Peto's method
#'     \item \code{"LR"}: Log-rank test based method
#'     \item \code{"Cox"}: Cox proportional hazards regression
#'   }
#' @param handling_tie A logical value indicating whether to handle tied observations.
#'   Set to FALSE for simulation data with continuous distributions (faster). Default is TRUE.
#' @param is_sorted A logical value indicating whether the data is already sorted by time.
#'   Set to TRUE for pre-sorted data to skip sorting operations (faster). Default is FALSE.
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{method}: The method used for estimation
#'     \item \code{HR}: The estimated hazard ratio
#'   }
#'
#' @details
#' This function implements five different methods for hazard ratio estimation:
#'
#' \describe{
#'   \item{Person-Year (PY)}{
#'     Estimates HR as the ratio of event rates:
#'     \deqn{HR = \frac{d_1/t_1}{d_0/t_0}}
#'     where \eqn{d_i} is the number of events and \eqn{t_i} is the total time at risk in group \eqn{i}.
#'   }
#'   \item{Pike's method}{
#'     Uses the formula: \deqn{HR = \frac{O_1 \times E_0}{O_0 \times E_1}}
#'     where \eqn{O_i} is observed events and \eqn{E_i} is expected events in group \eqn{i}.
#'   }
#'   \item{Peto's method}{
#'     Estimates HR as: \deqn{HR = \exp\left(\frac{O_1 - E_1}{V_1}\right)}
#'     where \eqn{V_1} is the variance of the log-rank statistic.
#'   }
#'   \item{Log-rank based (LR)}{
#'     Uses the log-rank statistic: \deqn{HR = \exp\left(\frac{(1+r)Z}{\sqrt{r(O_1 + O_0)}}\right)}
#'     where \eqn{r} is the allocation ratio and \eqn{Z} is the standardized log-rank statistic.
#'   }
#'   \item{Cox regression}{
#'     Uses the standard Cox proportional hazards model via \code{\link[survival]{coxph}}.
#'   }
#' }
#'
#' When \code{is_sorted = TRUE}, the function assumes the data is already sorted
#' by time and skips sorting operations for additional speed improvement.
#'
#' @examples
#' library(survival)
#'
#' # Compare different methods on ovarian cancer data
#' methods <- c('PY', 'Pike', 'Peto', 'LR', 'Cox')
#' results <- lapply(methods, function(m) {
#'   esthr(ovarian$futime, ovarian$fustat, ovarian$rx, 2, m)
#' })
#' hr_comparison <- do.call(rbind, results)
#' print(hr_comparison)
#'
#' # Using veteran lung cancer data (more events)
#' veteran_results <- lapply(methods, function(m) {
#'   esthr(veteran$time, veteran$status, veteran$trt, 1, m)
#' })
#' veteran_hr <- do.call(rbind, veteran_results)
#' print(veteran_hr)
#'
#' # Ultra-fast version for pre-sorted simulation data
#' set.seed(123)
#' n <- 1000
#' group_sim <- rep(c(1, 2), each = n/2)
#' time_sim <- c(rexp(n/2, 0.1), rexp(n/2, 0.15))
#' event_sim <- rep(1, n)
#' ord <- order(time_sim)
#' hr_fast <- esthr(time_sim[ord], event_sim[ord], group_sim[ord], 1, "Pike", is_sorted = TRUE)
#' print(hr_fast)
#'
#' @seealso
#' \code{\link[survival]{coxph}} for standard Cox regression,
#' \code{\link{lrtest}} for log-rank test statistics
#'
#' @importFrom survival coxph Surv
#' @importFrom stats coef
#' @export
esthr <- function(time, event, group, control, method, handling_tie = TRUE, is_sorted = FALSE) {
  # Input validation
  if (length(time) != length(event) || length(time) != length(group)) {
    stop("Arguments 'time', 'event', and 'group' must have the same length")
  }

  valid_methods <- c('PY', 'Pike', 'Peto', 'LR', 'Cox')
  if (!method %in% valid_methods) {
    stop(paste("Method must be one of:", paste(valid_methods, collapse = ", ")))
  }

  if (sum(event) == 0) {
    stop("No events observed in the data")
  }

  # Convert groups to binary indicator (0 = control, 1 = treatment)
  j <- as.numeric(group != control)

  # Method-specific calculations
  if (method == 'PY') {
    # Person-Year method - direct calculation without ranking
    rate_treatment <- sum(event * j) / sum(time * j)
    rate_control <- sum(event * (1 - j)) / sum(time * (1 - j))
    HR <- rate_treatment / rate_control

  } else if (method == 'Cox') {
    # Cox proportional hazards regression
    cox.fit <- coxph(Surv(time, event) ~ j)
    HR <- exp(coef(cox.fit)[1])

  } else {
    # Common calculations for Pike, Peto, and LR methods
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

      # Get event indices
      event_indices <- which(event_sorted == 1)
      n_events <- length(event_indices)
      n <- length(time)

      if (n_events == 0) {
        stop("No events observed in the data")
      }

      # Calculate numbers at risk at event times using efficient method
      n1k <- rev(cumsum(rev(j_sorted)))[event_indices]
      n0k <- (n:1)[event_indices] - n1k
      njk <- n1k + n0k

      # Events are single events at each time (no ties)
      e1k <- j_sorted[event_indices]
      e0k <- (1 - j_sorted)[event_indices]
      ejk <- rep(1, n_events)

      # Calculate observed and expected events
      O1 <- sum(e1k)
      O0 <- sum(e0k)
      E1 <- sum(ejk * (n1k / njk))
      E0 <- sum(ejk * (n0k / njk))

      # Method-specific hazard ratio calculations
      if (method == 'Pike') {
        HR <- (O1 * E0) / (O0 * E1)

      } else if (method == 'Peto') {
        V1 <- sum((n1k * n0k * ejk * (njk - ejk)) / (njk^2 * (njk - 1)), na.rm = TRUE)
        HR <- exp((O1 - E1) / V1)

      } else if (method == 'LR') {
        # Log-rank test based method with allocation ratio
        r <- sum(1 - j) / sum(j)  # Control/Treatment ratio
        V1 <- sum((n1k * n0k * ejk * (njk - ejk)) / (njk^2 * (njk - 1)), na.rm = TRUE)
        LR <- (O1 - E1) / sqrt(V1)
        HR <- exp((1 + r) * LR / sqrt(r * (O1 + O0)))
      }

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

      # Calculate numbers at risk using vectorized operations
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
      O0 <- sum(e.0k)
      E1 <- sum(e.jk * (n.1k / n.jk))
      E0 <- sum(e.jk * (n.0k / n.jk))

      # Method-specific hazard ratio calculations
      if (method == 'Pike') {
        HR <- (O1 * E0) / (O0 * E1)

      } else if (method == 'Peto') {
        V1 <- sum((n.1k * n.0k * e.jk * (n.jk - e.jk)) / (n.jk^2 * (n.jk - 1)), na.rm = TRUE)
        HR <- exp((O1 - E1) / V1)

      } else if (method == 'LR') {
        # Log-rank test based method with allocation ratio
        r <- sum(1 - j) / sum(j)  # Control/Treatment ratio
        V1 <- sum((n.1k * n.0k * e.jk * (n.jk - e.jk)) / (n.jk^2 * (n.jk - 1)), na.rm = TRUE)
        LR <- (O1 - E1) / sqrt(V1)
        HR <- exp((1 + r) * LR / sqrt(r * (O1 + O0)))
      }
    }
  }

  # Return result as data.frame
  result <- data.frame(
    method = method,
    HR = HR,
    stringsAsFactors = FALSE
  )
  return(result)
}
