#' Generate random numbers from piecewise exponential distribution
#'
#' This function generates random numbers from a piecewise exponential distribution
#' where each time interval has a different hazard rate. This implementation follows
#' the approach described in Luo et al. (2019) and is equivalent to the rpwe function
#' in the PWEALL package.
#'
#' @param n A positive integer specifying the number of observations to generate.
#' @param time A numeric vector of time points defining the intervals.
#'   The length of time should be one more than the length of hazard.
#'   The last element must be Inf.
#' @param hazard A numeric vector of hazard rates for each interval.
#'   The length should be one less than the length of time.
#'
#' @return A numeric vector of random numbers from the piecewise exponential distribution.
#'
#' @examples
#' # Standard exponential distribution (hazard = 1)
#' time1 <- c(0, Inf)
#' hazard1 <- 1
#' samples1 <- rpieceexp(1000, time1, hazard1)
#' hist(samples1, main = "Standard Exponential Distribution")
#'
#' # Piecewise exponential with multiple intervals
#' time2 <- c(0, 5, 10, Inf)
#' hazard2 <- c(0.1, 0.3, 0.2)
#' samples2 <- rpieceexp(1000, time2, hazard2)
#' hist(samples2, main = "Piecewise Exponential Distribution")
#'
#' # Three-piece exponential distribution
#' time3 <- c(0, 10, 20, Inf)
#' hazard3 <- c(0.2, 0.5, 0.1)
#' samples3 <- rpieceexp(1000, time3, hazard3)
#' hist(samples3, main = "Three-piece Exponential Distribution")
#'
#' @references
#' Luo, X., Mao, X., Chen, X., Qiu, J., Bai, S., & Quan, H. (2019).
#' Design and monitoring of survival trials in complex scenarios.
#' Statistics in Medicine, 38(2), 192-209.
#'
#' @import stats
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
