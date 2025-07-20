#' Simulate survival data for clinical trials using data.table
#'
#' This function generates simulation datasets for clinical trials with time-to-event endpoints
#' using data.table for efficient data manipulation while maintaining the high-performance
#' base R functions for random number generation. It accounts for piecewise uniform distribution
#' for patient accrual and piecewise exponential distributions for patient survival time and
#' dropout time during the trial.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 1000.
#' @param N A positive integer specifying the total sample size per simulation.
#' @param a.time A numeric vector of time points defining the accrual intervals.
#'   The length must be one more than the length of intensity (when used) or proportion (when used).
#' @param intensity A numeric vector of accrual intensities (rate per time unit) for each interval.
#'   The length should be one less than the length of a.time. Cannot be used with proportion.
#'   Default is NULL.
#' @param proportion A numeric vector of accrual proportions for each interval that sum to 1.
#'   The length should be one less than the length of a.time. Cannot be used with intensity.
#'   Default is NULL.
#' @param e.time A numeric vector of time points defining the survival time intervals.
#'   The length must be one more than the length of e.hazard. The last element must be Inf.
#' @param e.hazard A numeric vector of hazard rates for survival time in each interval.
#'   The length should be one less than the length of e.time.
#' @param d.time A numeric vector of time points defining the dropout time intervals.
#'   The length must be one more than the length of d.hazard. The last element must be Inf.
#' @param d.hazard A numeric vector of hazard rates for dropout time in each interval.
#'   The length should be one less than the length of d.time.
#' @param seed A positive integer for random seed setting. Default is NULL.
#'
#' @return A data.table containing the simulated survival data with the following columns:
#' \describe{
#'   \item{simID}{Simulation iteration ID (1 to nsim)}
#'   \item{accrual.time}{Patient accrual time from study start}
#'   \item{surv.time}{Survival time from patient entry}
#'   \item{dropout.time}{Time to dropout from patient entry}
#'   \item{tte}{Time-to-event (minimum of survival time and dropout time)}
#'   \item{total}{Total time from study start (accrual time + time-to-event)}
#'   \item{dropout}{Dropout indicator (1 = dropout occurred, 0 = event occurred)}
#' }
#'
#' @examples
#' library(data.table)
#'
#' # Basic simulation with uniform accrual
#' data1 <- simData(
#'   nsim = 100,
#'   N = 200,
#'   a.time = c(0, 24),
#'   intensity = 200/24,
#'   e.time = c(0, Inf),
#'   e.hazard = log(2) / 12,
#'   d.time = c(0, Inf),
#'   d.hazard = -log(1 - 0.1) / 12,
#'   seed = 123
#' )
#'
#' # Simulation with piecewise uniform accrual and piecewise exponential survival
#' data2 <- simData(
#'   nsim = 500,
#'   N = 300,
#'   a.time = c(0, 6, 12, 18, 24),
#'   intensity = c(10, 20, 30, 40),
#'   e.time = c(0, 6, Inf),
#'   e.hazard = c(0.1, 0.05),
#'   d.time = c(0, Inf),
#'   d.hazard = 0.02,
#'   seed = 456
#' )
#'
#' # Simulation using proportions for accrual
#' data3 <- simData(
#'   nsim = 1000,
#'   N = 400,
#'   a.time = c(0, 12, 24),
#'   proportion = c(0.3, 0.7),
#'   e.time = c(0, Inf),
#'   e.hazard = 0.08,
#'   d.time = c(0, Inf),
#'   d.hazard = 0.01
#' )
#'
#' @import data.table
#' @export
simData <- function(nsim = 1e+3, N, a.time, intensity = NULL, proportion = NULL, e.time,
                    e.hazard, d.time, d.hazard, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Input validation
  if (is.null(intensity) && is.null(proportion)) {
    stop("Either intensity or proportion must be specified")
  }

  if (!is.null(intensity) && !is.null(proportion)) {
    stop("Cannot specify both intensity and proportion arguments")
  }

  # Pre-calculate total sample size for efficient memory allocation
  total_n <- N * nsim

  # Create the main data.table with simID using efficient rep.int
  dt <- data.table(
    simID = rep.int(seq_len(nsim), rep.int(N, nsim))
  )

  # Generate accrual times using high-performance base function
  dt[, accrual.time := rpieceunif(.N, a.time, intensity, proportion)]

  # Generate survival times using high-performance base function
  dt[, surv.time := rpieceexp(.N, e.time, e.hazard)]

  # Generate dropout times using high-performance base function
  dt[, dropout.time := rpieceexp(.N, d.time, d.hazard)]

  # Calculate derived variables using data.table's efficient operations
  # Use := for reference modification (no copying)
  dt[, `:=`(
    tte = pmin(surv.time, dropout.time),
    dropout = as.integer(dropout.time < surv.time)
  )]

  # Calculate total time using data.table's efficient column operations
  dt[, total := accrual.time + tte]

  # Set column order efficiently without copying
  setcolorder(dt, c("simID", "accrual.time", "surv.time", "dropout.time", "tte", "total", "dropout"))

  return(dt)
}
