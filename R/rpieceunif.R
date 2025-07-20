#' Generate random numbers from piecewise uniform distribution
#'
#' This function generates random numbers from a piecewise uniform distribution
#' where each time interval has a different intensity (rate per time unit) or
#' a specified proportion of observations.
#'
#' @param n A positive integer specifying the number of observations to generate.
#' @param time A numeric vector of time points defining the intervals.
#'   The length of time should be one more than the length of intensity or proportion.
#' @param intensity A numeric vector of intensities (rate per time unit) for each interval.
#'   The length should be one less than the length of time. Cannot be used with proportion.
#' @param proportion A numeric vector of proportions for each interval that sum to 1.
#'   The length should be one less than the length of time. Cannot be used with intensity.
#'   Default is NULL.
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
#' # Generate 100 random numbers with specified proportions
#' time <- c(0, 6, 12, 18, 24)
#' proportion <- c(0.1, 0.2, 0.3, 0.4)
#' samples2 <- rpieceunif(100, time, proportion = proportion)
#' hist(samples2, main = "Piecewise Uniform with Proportions")
#'
#' # Example with decimal intensity
#' time2 <- c(0, 6)
#' intensity2 <- 24.5
#' samples3 <- rpieceunif(50, time2, intensity2)
#' hist(samples3, main = "Single Interval with Decimal Intensity")
#'
#' @import stats
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
