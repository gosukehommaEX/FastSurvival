#' Fast simulation data generation (C++ version)
#'
#' @param nsim Number of simulation iterations
#' @param N Total sample size per simulation
#' @param a.time Accrual time points
#' @param intensity Accrual intensities (or NULL)
#' @param proportion Accrual proportions (or NULL)
#' @param e.time Event time points
#' @param e.hazard Event hazards
#' @param d.time Dropout time points
#' @param d.hazard Dropout hazards
#' @param seed Random seed
#' @return data.table with simulation results
#' @export
simDataFast <- function(nsim = 1000, N, a.time, intensity = NULL, proportion = NULL,
                        e.time, e.hazard, d.time, d.hazard, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Input validation
  if (is.null(intensity) && is.null(proportion)) {
    stop("Either intensity or proportion must be specified")
  }

  # Convert arguments to appropriate types
  nsim <- as.integer(nsim)
  N <- as.integer(N)
  a.time <- as.numeric(a.time)
  e.time <- as.numeric(e.time)
  e.hazard <- as.numeric(e.hazard)
  d.time <- as.numeric(d.time)
  d.hazard <- as.numeric(d.hazard)

  # Call C++ function with corrected parameter order
  if (!is.null(intensity)) {
    intensity <- as.numeric(intensity)
    result <- simData_cpp(nsim, N, a.time, e.time, e.hazard, d.time, d.hazard,
                          intensity = intensity, proportion = NULL)
  } else {
    proportion <- as.numeric(proportion)
    result <- simData_cpp(nsim, N, a.time, e.time, e.hazard, d.time, d.hazard,
                          intensity = NULL, proportion = proportion)
  }

  # Convert to data.table
  library(data.table)
  dt <- data.table(
    simID = result$simID,
    accrual.time = result$accrual.time,
    surv.time = result$surv.time,
    dropout.time = result$dropout.time,
    tte = result$tte,
    total = result$total,
    dropout = result$dropout
  )

  return(dt)
}
