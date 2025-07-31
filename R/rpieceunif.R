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
#'
#' @return A numeric vector of random numbers from the piecewise uniform distribution.
#'
#' @examples
#' # Generate 100 random numbers with time intervals [0,6), [6,12), [12,18), [18,24)
#' # and intensities that sum to 100 total expected observations
#' time <- c(0, 6, 12, 18, 24)
#' intensity <- c(100/(6*4), 100/(6*4), 100/(6*4), 100/(6*4))  # Equal intensity per interval
#' samples <- rpieceunif(100, time, intensity)
#' hist(samples, main = "Piecewise Uniform Distribution")
#'
#' # Generate 100 random numbers with specified proportions
#' time <- c(0, 6, 12, 18, 24)
#' proportion <- c(0.1, 0.2, 0.3, 0.4)
#' samples2 <- rpieceunif(100, time, proportion = proportion)
#' hist(samples2, main = "Piecewise Uniform with Proportions")
#'
#' # Example with single interval where intensity * interval length = n
#' time2 <- c(0, 10)
#' intensity2 <- 5  # 5 * 10 = 50
#' samples3 <- rpieceunif(50, time2, intensity2)
#' hist(samples3, main = "Single Interval with Matching Expected Count")
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

  # Define tolerance for checking parameters
  tolerance <- .Machine$double.eps^0.5

  # Calculate cumulative proportions based on input type
  if (!is.null(proportion)) {
    # Input validation for proportion
    if (length(time) != length(proportion) + 1) {
      stop("Length of time must be one more than length of proportion")
    }

    if (any(proportion < 0)) {
      stop("All proportions must be non-negative")
    }

    if (abs(sum(proportion) - 1) > tolerance) {
      stop("Proportions must sum to 1")
    }

    cumulative_props <- cumsum(proportion)
  } else {
    # Input validation for intensity
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
    total_expected <- sum(expected_counts)

    # Check if total expected count matches requested n
    if (abs(total_expected - n) > tolerance) {
      stop(paste("Total expected count from intensity*intervals (",
                 round(total_expected, 6),
                 ") does not match requested n (", n, ")", sep = ""))
    }

    # Calculate cumulative proportions for each interval
    cumulative_props <- cumsum(expected_counts) / total_expected
  }

  # Generate random numbers for interval selection
  u_interval <- runif(n)

  # Vectorized interval assignment
  interval_idx <- as.numeric(cut(u_interval, breaks = c(0, cumulative_props),
                                 labels = FALSE, include.lowest = TRUE))

  # Generate uniform random numbers within assigned intervals
  # This requires NEW independent uniform random numbers
  u_position <- runif(n)
  start_times <- time[interval_idx]
  end_times <- time[interval_idx + 1]
  result <- start_times + u_position * (end_times - start_times)

  return(result)
}
