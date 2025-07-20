#' Simulate survival data for clinical trials with time-to-event endpoints
#'
#' This function generates datasets for clinical trials with time-to-event endpoints.
#' The function can account for multiple subgroups, multiple analysis timing for
#' interim analyses, piecewise uniform distribution for patient accrual, and dropout
#' during the trial.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 10000.
#' @param mst.T A numeric value representing the median survival time for the treatment group.
#' @param mst.C A numeric value representing the median survival time for the control group.
#' @param N.T A numeric vector specifying the sample size of the treatment group in each subgroup.
#' @param N.C A numeric vector specifying the sample size of the control group in each subgroup.
#' @param a.time A numeric vector of time points defining the accrual intervals.
#' @param a.intensity A numeric vector of intensities (rate per time unit) for each accrual interval.
#' @param d.rate A numeric value representing the dropout rate per year. Default is 0.
#' @param E A numeric vector specifying the number of events for each interim analysis.
#' @param seed A positive integer for random seed setting. Default is 1.
#'
#' @return A tibble containing the simulated survival data with the following columns:
#' \describe{
#'   \item{sim}{Simulation iteration ID}
#'   \item{subgroup}{Subgroup ID}
#'   \item{trt}{Treatment indicator (1 = treatment, 0 = control)}
#'   \item{accrual}{Patient accrual time}
#'   \item{followup}{Follow-up time}
#'   \item{dropout}{Time at dropout}
#'   \item{tte}{Time-to-event (minimum of follow-up and dropout time)}
#'   \item{total}{Total time (accrual + time-to-event)}
#'   \item{dropout.flag}{Dropout indicator (1 = dropout, 0 = no dropout)}
#'   \item{event}{Event indicator matrix for interim analyses (1 = event, 0 = censored)}
#' }
#'
#' @examples
#' # Basic simulation with multiple subgroups and interim analyses
#' data <- simSurvdata(
#'   nsim = 100,
#'   mst.T = 5.811,
#'   mst.C = 4.3,
#'   N.T = c(25, 112, 113),
#'   N.C = c(25, 112, 113),
#'   a.time = c(0, 5, 10, 15, 20, 25),
#'   a.intensity = c(3.5, 14.3, 28.9, 43.6, 45),
#'   d.rate = 0.01,
#'   E = c(142, 248, 354),
#'   seed = 123
#' )
#'
#' # Simple two-arm trial simulation
#' data_simple <- simSurvdata(
#'   nsim = 50,
#'   mst.T = 12,
#'   mst.C = 8,
#'   N.T = 100,
#'   N.C = 100,
#'   a.time = c(0, 12),
#'   a.intensity = 20,
#'   d.rate = 0.05,
#'   E = c(80, 120),
#'   seed = 456
#' )
#'
#' @import dplyr
#' @import tidyr
#' @export
simSurvdata <- function(nsim = 1e+4, mst.T, mst.C, N.T, N.C, a.time, a.intensity, d.rate = 0, E, seed = 1) {

  # Set seed number
  set.seed(seed)

  # Number of subgroups
  if (length(N.T) != length(N.C)) {
    stop('Number of subgroups must be the same between treatment and control groups')
  }
  Subg <- length(N.T)

  # Sample size
  N.js <- rbind(N.T, N.C)
  N <- sum(N.js)

  # Generate dataset
  dataset <- tibble(
    sim = rep(1:nsim, rep(N, nsim)), # Simulation ID
    subgroup = rep(rep(rep(1:Subg, 2), as.vector(t(N.js))), nsim), # Subgroup ID
    trt = rep(rep(1:0, rowSums(N.js)), nsim), # Treatment ID (1 = treatment, 0 = control)
    accrual = rpieceunif(N * nsim, a.time, a.intensity), # Accrual time
    followup = (trt == 1) * rexp(N * nsim, rate = log(2) / mst.T) + (trt == 0) * rexp(N * nsim, rate = log(2) / mst.C), # Follow-up time
    dropout = rexp(N * nsim, rate = -log(1 - d.rate) / 12), # Time at dropout
    tte = pmin(followup, dropout), # Time-to-event
    total = accrual + tte, # Total time (accrual + time-to-event)
    dropout.flag = as.numeric(dropout < followup) # Dropout indicator
  ) %>%
    arrange(sim, total) %>%
    group_by(sim) %>%
    mutate(
      tte = '+'(
        outer(total, total[dropout.flag != 1][E], FUN = '>') * pmax(0, outer(-accrual, total[dropout.flag != 1][E], FUN = '+')),
        outer(total, total[dropout.flag != 1][E], FUN = '<=') * tte
      ),
      event = outer(total, total[dropout.flag != 1][E], FUN = '<=') * (dropout.flag != 1) # Event indicator for interim analyses
    ) %>%
    ungroup()

  # Return the dataset
  return(dataset)
}
