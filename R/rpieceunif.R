#' Generate random numbers from piecewise uniform distribution
#'
#' This function generates random numbers from a piecewise uniform distribution
#' where each time interval has a different intensity (rate per time unit) or
#' a specified proportion of observations. This is particularly useful for modeling
#' non-uniform recruitment patterns in clinical trials.
#'
#' @param n A positive integer specifying the number of observations to generate.
#' @param time A numeric vector of time points defining the intervals.
#'   The length of time should be one more than the length of intensity or proportion.
#'   Time points must be in increasing order.
#' @param intensity A numeric vector of intensities (rate per time unit) for each interval.
#'   The length should be one less than the length of time. Cannot be used with proportion.
#'   All values must be positive. Default is NULL.
#' @param proportion A numeric vector of proportions for each interval that sum to 1.
#'   The length should be one less than the length of time. Cannot be used with intensity.
#'   All values must be non-negative and sum to 1. Default is NULL.
#'
#' @return A numeric vector of random numbers from the piecewise uniform distribution.
#'
#' @details
#' The piecewise uniform distribution allows for non-constant recruitment rates
#' across different time intervals. This is commonly used in clinical trial
#' simulations where patient recruitment patterns vary over time.
#'
#' When using \code{intensity}, the probability of recruitment in interval \eqn{i}
#' is proportional to \eqn{intensity_i \times (time_{i+1} - time_i)}.
#'
#' When using \code{proportion}, the values directly specify the probability
#' of recruitment in each interval.
#'
#' @examples
#' # Generate 1000 random numbers with varying intensities
#' time <- c(0, 6, 12, 18, 24)
#' intensity <- c(10, 20, 30, 40)
#' samples <- rpieceunif(1000, time, intensity = intensity)
#' hist(samples, main = "Piecewise Uniform Distribution (Intensity)",
#'      xlab = "Time", breaks = 20)
#'
#' # Generate 1000 random numbers with specified proportions
#' time <- c(0, 6, 12, 18, 24)
#' proportion <- c(0.1, 0.2, 0.3, 0.4)
#' samples2 <- rpieceunif(1000, time, proportion = proportion)
#' hist(samples2, main = "Piecewise Uniform Distribution (Proportions)",
#'      xlab = "Time", breaks = 20)
#'
#' # Simple example with single interval
#' time_simple <- c(0, 10)
#' intensity_simple <- 5
#' samples3 <- rpieceunif(500, time_simple, intensity = intensity_simple)
#' hist(samples3, main = "Single Interval Uniform Distribution",
#'      xlab = "Time", breaks = 15)
#'
#' # Clinical trial recruitment example
#' # Slow start, ramp up, then slow down
#' trial_time <- c(0, 3, 9, 15, 18)
#' trial_intensity <- c(2, 8, 12, 4)  # patients per month
#' recruitment <- rpieceunif(200, trial_time, intensity = trial_intensity)
#' plot(sort(recruitment), 1:200, type = "s",
#'      xlab = "Time (months)", ylab = "Cumulative Patients",
#'      main = "Clinical Trial Recruitment Pattern")
#'
#' @seealso
#' \code{\link{rpieceexp}} for piecewise exponential distribution,
#' \code{\link{simData}} for survival data simulation using piecewise distributions
#'
#' @references
#' Luo, X., Mao, X., Chen, X., Qiu, J., Bai, S., & Quan, H. (2019).
#' Design and monitoring of survival trials in complex scenarios.
#' Statistics in Medicine, 38(2), 192-209.
#'
#' @importFrom stats runif
#' @export
rpieceunif <- function(n, time, intensity = NULL, proportion = NULL) {
  # Input validation for mutual exclusivity
  if (is.null(intensity) && is.null(proportion)) {
    stop("Either intensity or proportion must be specified")
  }

  if (!is.null(intensity) && !is.null(proportion)) {
    stop("Cannot specify both intensity and proportion arguments")
  }

  # Input validation for proportion
  if (!is.null(proportion)) {
    if (length(time) != length(proportion) + 1) {
      stop("Length of time must be one more than length of proportion")
    }

    if (any(proportion < 0)) {
      stop("All proportions must be non-negative")
    }

    if (abs(sum(proportion) - 1) > .Machine$double.eps^0.5) {
      stop("Proportions must sum to 1")
    }

    # Calculate cumulative proportions for each interval
    cumulative_props <- cumsum(proportion)

    # Generate random numbers for interval selection
    u <- runif(n)

    # Vectorized interval assignment
    # Use cut() to assign intervals based on cumulative proportions
    interval_idx <- as.numeric(cut(u, breaks = c(0, cumulative_props), labels = FALSE, include.lowest = TRUE))

    # Generate uniform random numbers within assigned intervals (vectorized)
    start_times <- time[interval_idx]
    end_times <- time[interval_idx + 1]
    result <- runif(n, min = start_times, max = end_times)

    return(result)
  }

  # Original intensity-based implementation
  if (length(time) != length(intensity) + 1) {
    stop("Length of time must be one more than length of intensity")
  }

  if (any(intensity <= 0)) {
    stop("All intensities must be positive")
  }

  if (any(diff(time) <= 0)) {
    stop("Time vector must be strictly increasing")
  }

  # Calculate interval lengths
  interval_lengths <- diff(time)

  # Calculate total expected number of observations in each interval
  expected_counts <- intensity * interval_lengths

  # Calculate cumulative proportions for each interval
  total_expected <- sum(expected_counts)
  cumulative_props <- cumsum(expected_counts) / total_expected

  # Generate random numbers for interval selection
  u <- runif(n)

  # Vectorized interval assignment
  # Use cut() to assign intervals based on cumulative proportions
  interval_idx <- as.numeric(cut(u, breaks = c(0, cumulative_props), labels = FALSE, include.lowest = TRUE))

  # Generate uniform random numbers within assigned intervals (vectorized)
  start_times <- time[interval_idx]
  end_times <- time[interval_idx + 1]
  result <- runif(n, min = start_times, max = end_times)

  return(result)
}
