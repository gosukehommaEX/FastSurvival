#' Simulate survival data for clinical trials with time-to-event endpoints
#'
#' This function generates datasets for clinical trials with time-to-event endpoints.
#' The function can account for piecewise uniform distribution for patient accrual, and piecewise
#' exponential distributions for patient survival time and dropout during the trial.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 10000.
#' @param N A numeric vector specifying the sample size.
#' @param a.time A numeric vector of time points defining the accrual intervals.
#' @param intensity A numeric vector of intensities (rate per time unit) for each interval.
#'   The length should be one less than the length of time. Cannot be used with proportion.
#' @param proportion A numeric vector of proportions for each interval that sum to 1.
#'   The length should be one less than the length of time. Cannot be used with intensity.
#'   Default is NULL.
#' @param e.time A numeric vector of time points defining the intervals.
#'   The length of time should be one more than the length of hazard.
#'   The last element must be Inf.
#' @param e.hazard A numeric vector of hazard rates for each interval.
#'   The length should be one less than the length of time.
#' @param d.time A numeric vector of time points defining the intervals.
#'   The length of time should be one more than the length of hazard.
#'   The last element must be Inf.
#' @param d.hazard A numeric vector of hazard rates for each interval.
#'   The length should be one less than the length of time.
#' @param seed A positive integer for random seed setting. Default is NULL.
#'
#' @return A tibble containing the simulated survival data with the following columns:
#' \describe{
#'   \item{simID}{Simulation iteration ID}
#'   \item{accrual.time}{Patient accrual time}
#'   \item{surv.time}{Survival time}
#'   \item{dropout.time}{Time at dropout}
#'   \item{tte}{Time-to-event (minimum of follow-up and dropout time)}
#'   \item{total}{Total time (accrual + time-to-event)}
#'   \item{dropout}{Dropout indicator (1 = dropout, 0 = no dropout)}
#' }
#'
#' @examples
#' # Basic simulation
#' data <- simSurvdata(
#'   nsim = 1e+3,
#'   N = 200,
#'   a.time = c(0, 6, 12, 18, 24),
#'   intensity = c(10, 20, 30, 40),
#'   proportion = NULL,
#'   e.time = c(0, Inf),
#'   e.hazard = log(2) / 12,
#'   d.time = c(0, Inf),
#'   d.hazard = -log(1 - 0.1) / 12,
#'   seed = NULL
#' )
#'
#' @import dplyr
#' @export
simSurvdata <- function(nsim = 1e+3, N, a.time, intensity, proportion = NULL, e.time,
                        e.hazard, d.time, d.hazard, seed = NULL) {

  # Set seed number
  set.seed(seed)

  # Generate dataset
  dataset <- tibble(
    simID = rep(1:nsim, rep(N, nsim)), # Simulation ID
    accrual.time = rpieceunif(N * nsim, a.time, intensity, proportion), # Accrual time
    surv.time = rpieceexp(N * nsim, e.time, e.hazard), # Follow-up time
    dropout.time = rpieceexp(N * nsim, d.time, d.hazard), # Time at dropout
    tte = pmin(surv.time, dropout.time), # Time-to-event
    total = accrual.time + tte, # Total time (accrual + time-to-event)
    dropout = as.numeric(dropout.time < surv.time) # Dropout indicator
  )

  # Return the dataset
  return(dataset)
}
