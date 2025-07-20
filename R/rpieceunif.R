#' Generate random numbers from piecewise uniform distribution
#'
#' This function generates random numbers from a piecewise uniform distribution
#' where each time interval has a different intensity (rate per time unit).
#'
#' @param n A positive integer specifying the number of observations to generate.
#' @param time A numeric vector of time points defining the intervals.
#'   The length of time should be one more than the length of intensity.
#' @param intensity A numeric vector of intensities (rate per time unit) for each interval.
#'   The length should be one less than the length of time.
#'
#' @return A numeric vector of random numbers from the piecewise uniform distribution.
#'
#' @examples
#' # Generate 100 random numbers with time intervals [0,6), [6,12), [12,18), [18,24)
#' # and intensities 10, 20, 30, 40 per time unit respectively
#' time <- c(0, 6, 12, 18, 24)
#' intensity <- c(10, 20, 30, 40)
#' samples <- rpieceunif(100, time, intensity)
#' hist(samples, main = "Piecewise Uniform Distribution")
#'
#' # Example with decimal intensity
#' time2 <- c(0, 6)
#' intensity2 <- 24.5
#' samples2 <- rpieceunif(50, time2, intensity2)
#' hist(samples2, main = "Single Interval with Decimal Intensity")
#'
#' @import stats
#' @export
rpieceunif <- function(n, time, intensity) {
  # Input validation
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
