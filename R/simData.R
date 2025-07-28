#' Simulate survival data for clinical trials using tidyverse
#'
#' This function generates simulation datasets for clinical trials with time-to-event
#' endpoints using tidyverse packages for efficient data manipulation. It accounts for
#' piecewise uniform distribution for patient accrual and piecewise exponential
#' distributions for patient survival time and dropout time during the trial.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 1000.
#' @param n A positive integer specifying the total sample size per simulation.
#' @param a.time A numeric vector of time points defining the accrual intervals.
#'   The length must be one more than the length of intensity (when used) or proportion (when used).
#'   Time points must be in increasing order.
#' @param intensity A numeric vector of accrual intensities (rate per time unit) for each interval.
#'   The length should be one less than the length of a.time. Cannot be used with proportion.
#'   All values must be positive. Default is NULL.
#' @param proportion A numeric vector of accrual proportions for each interval that sum to 1.
#'   The length should be one less than the length of a.time. Cannot be used with intensity.
#'   All values must be non-negative and sum to 1. Default is NULL.
#' @param e.time A numeric vector of time points defining the survival time intervals.
#'   The length must be one more than the length of e.hazard. The last element must be Inf.
#'   Time points must be in increasing order.
#' @param e.hazard A numeric vector of hazard rates for survival time in each interval.
#'   The length should be one less than the length of e.time. All values must be positive.
#' @param d.time A numeric vector of time points defining the dropout time intervals.
#'   The length must be one more than the length of d.hazard. The last element must be Inf.
#'   Time points must be in increasing order. Can be NULL if no dropout is assumed. Default is NULL.
#' @param d.hazard A numeric vector of hazard rates for dropout time in each interval.
#'   The length should be one less than the length of d.time. All values must be positive.
#'   Can be NULL if no dropout is assumed. Default is NULL.
#' @param seed A positive integer for random seed setting. Default is NULL.
#'
#' @return A tibble containing the simulated survival data with the following columns:
#' \describe{
#'   \item{simID}{Simulation iteration ID (1 to nsim)}
#'   \item{patientID}{Patient ID (1 to n)}
#'   \item{accrual.time}{Patient accrual time from study start}
#'   \item{surv.time}{Survival time from patient entry}
#'   \item{dropout.time}{Time to dropout from patient entry}
#'   \item{tte}{Time-to-event (minimum of survival time and dropout time)}
#'   \item{total}{Total time from study start (accrual time + time-to-event)}
#'   \item{dropout}{Dropout indicator (1 = dropout occurred, 0 = event occurred)}
#' }
#'
#' @details
#' This function is designed for single-arm or single-group survival data simulation.
#' For multi-arm trials, use \code{\link{simTrial}}. The function uses piecewise
#' distributions to model:
#'
#' \describe{
#'   \item{Patient Accrual}{
#'     Modeled using piecewise uniform distribution with either intensities or proportions.
#'     This allows for realistic recruitment patterns that vary over time.
#'   }
#'   \item{Survival Times}{
#'     Modeled using piecewise exponential distribution to accommodate scenarios such as
#'     delayed treatment effects or changing baseline risks over time.
#'   }
#'   \item{Dropout Times}{
#'     Also modeled using piecewise exponential distribution to account for
#'     time-varying dropout rates. Can be disabled by setting d.time and d.hazard to NULL.
#'   }
#' }
#'
#' The observed time-to-event is the minimum of survival time and dropout time.
#' The total time from study start combines accrual time and time-to-event.
#' When dropout is not assumed, dropout.time is set to Inf and dropout indicator is 0.
#'
#' @examples
#' library(dplyr)
#'
#' # Basic simulation with uniform accrual over 24 months
#' data1 <- simData(
#'   nsim = 100,
#'   n = 200,
#'   a.time = c(0, 24),
#'   intensity = 200/24,
#'   e.time = c(0, Inf),
#'   e.hazard = log(2) / 12,
#'   d.time = c(0, Inf),
#'   d.hazard = -log(1 - 0.1) / 12,
#'   seed = 123
#' )
#' print(data1)
#'
#' # Simulation without dropout
#' data2 <- simData(
#'   nsim = 50,
#'   n = 300,
#'   a.time = c(0, 6, 12, 18, 24),
#'   intensity = c(5, 15, 25, 10),
#'   e.time = c(0, 6, Inf),
#'   e.hazard = c(0.1, 0.05),
#'   d.time = NULL,
#'   d.hazard = NULL,
#'   seed = 456
#' )
#'
#' @seealso
#' \code{\link{simTrial}} for multi-arm clinical trial simulation,
#' \code{\link{rpieceexp}} and \code{\link{rpieceunif}} for the underlying distributions,
#' \code{\link{analysisData}} for creating analysis datasets
#'
#' @import dplyr
#' @importFrom stats runif
#' @export
simData <- function(nsim = 1000, n, a.time, intensity = NULL, proportion = NULL,
                    e.time, e.hazard, d.time = NULL, d.hazard = NULL, seed = NULL) {

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

  # Validate dropout parameters
  if (is.null(d.time) && !is.null(d.hazard)) {
    stop("If d.time is NULL, d.hazard must also be NULL")
  }

  if (!is.null(d.time) && is.null(d.hazard)) {
    stop("If d.hazard is NULL, d.time must also be NULL")
  }

  # Check if dropout is assumed
  has_dropout <- !is.null(d.time) && !is.null(d.hazard)

  # Create the main dataset using tidyverse approach
  dataset <- tibble(
    simID = rep(1:nsim, each = n),
    patientID = rep(1:n, nsim)
  ) %>%
    group_by(simID) %>%
    mutate(
      # Generate accrual times using rpieceunif
      accrual.time = rpieceunif(n, a.time, intensity, proportion),
      # Generate survival times using rpieceexp
      surv.time = rpieceexp(n, e.time, e.hazard),
      # Generate dropout times based on whether dropout is assumed
      dropout.time = if (has_dropout) {
        rpieceexp(n, d.time, d.hazard)
      } else {
        rep(Inf, n)
      }
    ) %>%
    mutate(
      # Calculate time-to-event and total time
      tte = pmin(surv.time, dropout.time),
      total = accrual.time + tte,
      dropout = as.numeric(dropout.time < surv.time)
    ) %>%
    ungroup()

  return(dataset)
}
