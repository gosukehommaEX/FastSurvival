#' Generate random numbers from piecewise Poisson distribution
#'
#' This function generates random numbers from a piecewise Poisson distribution
#' where each time interval has a different rate parameter. This is particularly
#' useful for modeling event arrivals or patient enrollment in clinical trials
#' where the arrival rate changes over time.
#'
#' @param n A positive integer specifying the number of observations to generate.
#' @param time A numeric vector of time points defining the start of each interval.
#'   The length of time should be equal to the length of rate.
#'   Time points must be in strictly increasing order. The final interval
#'   extends to infinity with the last rate.
#' @param rate A numeric vector of rate parameters (events per time unit) for each interval.
#'   The length should be equal to the length of time.
#'   All rates must be non-negative.
#' @param continuous A logical value indicating the type of time generation.
#'   If TRUE, generates continuous event times within each time unit using uniform
#'   distribution (realistic for arrival processes). If FALSE, generates discrete
#'   integer event times (appropriate for daily/monthly counting). Default is FALSE.
#'
#' @return A numeric vector of random numbers from the piecewise Poisson distribution.
#'   If continuous = TRUE, returns real-valued event times sorted in ascending order.
#'   If continuous = FALSE, returns integer-valued event times representing discrete
#'   time points where events occur.
#'
#' @details
#' The piecewise Poisson distribution allows for modeling event arrival processes
#' where the rate changes at pre-specified time points. This is commonly used in
#' clinical trial simulations for patient enrollment modeling.
#'
#' The rate function is defined as:
#' \deqn{\lambda(t) = \lambda_j \quad \text{for } t_{j-1} \leq t < t_j}
#' where \eqn{\lambda_j} is the rate parameter in interval \eqn{j}.
#'
#' **Discrete Mode (continuous = FALSE)**:
#' For each integer time point t in interval j, generates Poisson(λ_j) events
#' occurring exactly at time t. This mimics the original discrete-time approach
#' where events are counted at regular intervals.
#'
#' **Continuous Mode (continuous = TRUE)**:
#' For each time unit [t, t+1) in interval j, generates Poisson(λ_j) events
#' and places them uniformly within the time unit. This provides a more realistic
#' continuous-time Poisson process approximation.
#'
#' The algorithm uses an efficient vectorized approach that:
#' \itemize{
#'   \item Calculates expected counts per time unit for each interval
#'   \item Generates Poisson counts for all time units simultaneously
#'   \item Handles discrete or continuous positioning based on the continuous parameter
#'   \item Extends generation beyond specified intervals if needed to reach sample size n
#' }
#'
#' @examples
#' # Generate 100 arrivals with constant rate (discrete)
#' time1 <- c(0)
#' rate1 <- c(5)
#' samples1_discrete <- rpiecepois(100, time1, rate1, continuous = FALSE)
#' hist(samples1_discrete, main = "Discrete Poisson Process", xlab = "Time", breaks = 20)
#'
#' # Generate 100 arrivals with constant rate (continuous)
#' samples1_continuous <- rpiecepois(100, time1, rate1, continuous = TRUE)
#' hist(samples1_continuous, main = "Continuous Poisson Process", xlab = "Time", breaks = 20)
#'
#' # Generate 650 arrivals with varying rates over different periods
#' time2 <- c(0, 5, 10, 15, 20)
#' rate2 <- c(3.5, 14.3, 28.9, 43.6, 45)
#' samples2 <- rpiecepois(650, time2, rate2, continuous = TRUE)
#' plot(samples2, seq_along(samples2), type = "s",
#'      xlab = "Time", ylab = "Cumulative Count",
#'      main = "Piecewise Poisson Enrollment Pattern")
#'
#' # Clinical trial enrollment with discrete daily counts
#' time_trial <- c(0, 3, 12, 18, 24)
#' rate_trial <- c(2, 10, 15, 8, 5)  # Events per day, final rate = 5 for t >= 24
#' enrollment_discrete <- rpiecepois(300, time_trial, rate_trial, continuous = FALSE)
#' table(enrollment_discrete)  # Count events per day
#'
#' # Clinical trial enrollment with continuous arrival times
#' enrollment_continuous <- rpiecepois(300, time_trial, rate_trial, continuous = TRUE)
#' hist(enrollment_continuous, main = "Continuous Enrollment Times",
#'      xlab = "Time (days)", breaks = 30)
#'
#' @importFrom stats rpois runif
#' @export
rpiecepois <- function(n, time, rate, continuous = FALSE) {
  # Input validation
  if (length(time) != length(rate)) {
    stop("Length of time must be equal to length of rate")
  }

  if (any(rate < 0)) {
    stop("All rates must be non-negative")
  }

  if (any(diff(time) <= 0)) {
    stop("Time vector must be strictly increasing")
  }

  if (n <= 0) {
    stop("n must be a positive integer")
  }

  # Create intervals: time[i] to time[i+1] (or infinity for last interval)
  n_intervals <- length(time)

  # Vectorized approach: generate all time units and their Poisson counts
  all_times <- numeric(0)
  all_counts <- numeric(0)

  # Process defined intervals
  for (i in seq_len(n_intervals - 1)) {
    period_start <- time[i]
    period_end <- time[i + 1]
    period_rate <- rate[i]

    if (period_rate > 0) {
      # Create time sequence for this period
      time_units <- seq(from = period_start, to = period_end - 1, by = 1)

      # Generate Poisson counts for all time units in this period
      poisson_counts <- rpois(length(time_units), period_rate)

      # Store results
      all_times <- c(all_times, time_units)
      all_counts <- c(all_counts, poisson_counts)
    }
  }

  # Generate individual event times based on counts and continuous parameter
  event_times <- numeric(0)

  for (i in seq_along(all_times)) {
    if (all_counts[i] > 0) {
      if (continuous) {
        # Continuous mode: uniform random times within [t, t+1)
        unit_times <- all_times[i] + runif(all_counts[i])
        event_times <- c(event_times, unit_times)
      } else {
        # Discrete mode: exact integer times
        unit_times <- rep(all_times[i], all_counts[i])
        event_times <- c(event_times, unit_times)
      }
    }
  }

  # Handle case where we don't have enough events using the last rate
  if (length(event_times) < n) {
    last_rate <- rate[n_intervals]  # Use the last rate for extension
    current_time <- time[n_intervals]  # Start from the last time point

    while (length(event_times) < n) {
      additional_count <- rpois(1, last_rate)
      if (additional_count > 0) {
        if (continuous) {
          additional_times <- current_time + runif(additional_count)
        } else {
          additional_times <- rep(current_time, additional_count)
        }
        event_times <- c(event_times, additional_times)
      }
      current_time <- current_time + 1
    }
  }

  # Sort and truncate to exact sample size
  event_times <- sort(event_times)[1:n]

  return(event_times)
}
