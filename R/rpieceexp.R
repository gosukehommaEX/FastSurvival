#' Generate random numbers from piecewise exponential distribution
#'
#' This function generates random numbers from a piecewise exponential distribution
#' where each time interval has a different hazard rate. This implementation follows
#' the approach described in Luo et al. (2019) and is commonly used for modeling
#' time-to-event data with non-constant hazard rates over time.
#'
#' @param n A positive integer specifying the number of observations to generate.
#' @param time A numeric vector of time points defining the intervals.
#'   The length of time should be one more than the length of hazard.
#'   The last element must be Inf to represent the final open interval.
#'   Time points must be in strictly increasing order.
#' @param hazard A numeric vector of hazard rates for each interval.
#'   The length should be one less than the length of time.
#'   All hazard rates must be positive.
#'
#' @return A numeric vector of random numbers from the piecewise exponential distribution.
#'
#' @details
#' The piecewise exponential distribution allows for modeling survival times where
#' the hazard rate changes at pre-specified time points. This is particularly useful
#' for clinical trials where treatment effects may be delayed or where the underlying
#' risk changes over time.
#'
#' The hazard function is defined as:
#' \deqn{h(t) = h_j \quad \text{for } t_{j-1} \leq t < t_j}
#' where \eqn{h_j} is the hazard rate in interval \eqn{j}.
#'
#' The survival function is:
#' \deqn{S(t) = \exp\left(-\sum_{j=1}^{k-1} h_j (t_j - t_{j-1}) - h_k (t - t_{k-1})\right)}
#' where \eqn{t_{k-1} \leq t < t_k}.
#'
#' @examples
#' # Standard exponential distribution (single interval)
#' time1 <- c(0, Inf)
#' hazard1 <- 0.1
#' samples1 <- rpieceexp(1000, time1, hazard1)
#' hist(samples1, main = "Standard Exponential Distribution",
#'      xlab = "Time", probability = TRUE, breaks = 30)
#' # Overlay theoretical density
#' x <- seq(0, max(samples1), length.out = 100)
#' lines(x, 0.1 * exp(-0.1 * x), col = "red", lwd = 2)
#'
#' # Piecewise exponential with increasing hazard
#' time2 <- c(0, 5, 10, Inf)
#' hazard2 <- c(0.05, 0.15, 0.25)  # Increasing hazard over time
#' samples2 <- rpieceexp(1000, time2, hazard2)
#' hist(samples2, main = "Piecewise Exponential (Increasing Hazard)",
#'      xlab = "Time", probability = TRUE, breaks = 30)
#'
#' # Piecewise exponential with decreasing hazard (improving prognosis)
#' time3 <- c(0, 2, 6, 12, Inf)
#' hazard3 <- c(0.8, 0.4, 0.2, 0.1)  # Decreasing hazard
#' samples3 <- rpieceexp(1000, time3, hazard3)
#' hist(samples3, main = "Piecewise Exponential (Decreasing Hazard)",
#'      xlab = "Time", probability = TRUE, breaks = 30)
#'
#' # Clinical trial example: delayed treatment effect
#' # High initial hazard, then lower after treatment kicks in
#' time_trial <- c(0, 3, Inf)  # Treatment effect starts at 3 months
#' hazard_control <- c(0.2, 0.2)  # Constant hazard for control
#' hazard_treatment <- c(0.2, 0.1)  # Reduced hazard after 3 months
#'
#' control_times <- rpieceexp(500, time_trial, hazard_control)
#' treatment_times <- rpieceexp(500, time_trial, hazard_treatment)
#'
#' # Compare survival curves
#' library(survival)
#' combined_data <- data.frame(
#'   time = c(control_times, treatment_times),
#'   event = rep(1, 1000),
#'   group = rep(c("Control", "Treatment"), each = 500)
#' )
#'
#' plot(survfit(Surv(time, event) ~ group, data = combined_data),
#'      main = "Survival Curves with Delayed Treatment Effect",
#'      xlab = "Time", ylab = "Survival Probability",
#'      col = c("red", "blue"))
#' legend("topright", c("Control", "Treatment"), col = c("red", "blue"), lty = 1)
#'
#' @seealso
#' \code{\link{rpieceunif}} for piecewise uniform distribution,
#' \code{\link{simData}} for survival data simulation,
#' \code{\link[stats]{rexp}} for standard exponential distribution
#'
#' @references
#' Luo, X., Mao, X., Chen, X., Qiu, J., Bai, S., & Quan, H. (2019).
#' Design and monitoring of survival trials in complex scenarios.
#' Statistics in Medicine, 38(2), 192-209.
#'
#' @importFrom stats runif
#' @export
rpieceexp <- function(n, time, hazard) {
  # Input validation
  if (length(time) != length(hazard) + 1) {
    stop("Length of time must be one more than length of hazard")
  }

  if (any(hazard <= 0)) {
    stop("All hazard rates must be positive")
  }

  if (!is.infinite(time[length(time)])) {
    stop("Last element of time must be Inf")
  }

  if (any(diff(time[1:(length(time)-1)]) <= 0)) {
    stop("Finite time points must be strictly increasing")
  }

  # Generate uniform random numbers
  u <- runif(n)

  # Calculate cumulative hazard at time change points
  # For each finite interval i: length = time[i+1] - time[i]
  n_intervals <- length(hazard)
  interval_lengths <- diff(time[1:n_intervals])  # Only finite intervals

  # Cumulative hazard: H(t_i) = sum_{j=1}^{i-1} hazard_j * length_j
  cumulative_hazard <- c(0, cumsum(hazard[1:(n_intervals-1)] * interval_lengths))

  # Use inverse CDF method
  target_hazard <- -log(u)

  # Find which interval each target falls into
  interval_idx <- findInterval(target_hazard, cumulative_hazard, rightmost.closed = TRUE)
  interval_idx <- pmax(1, interval_idx)
  interval_idx <- pmin(interval_idx, n_intervals)

  # Vectorized calculation
  start_times <- time[interval_idx]
  current_hazard <- hazard[interval_idx]
  prev_cumulative_hazard <- cumulative_hazard[interval_idx]

  remaining_hazard <- target_hazard - prev_cumulative_hazard
  result <- start_times + remaining_hazard / current_hazard

  return(result)
}
