#' Simulate clinical trial data for multiple groups
#'
#' This function generates simulation datasets for clinical trials with multiple treatment
#' groups using group-wise operations. It supports flexible accrual patterns (piecewise
#' uniform or piecewise Poisson) and optional dropout modeling.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 1000.
#' @param n A named list where each element corresponds to a group with sample sizes.
#'   Example: list(control = 100, treatment = 100) for two groups.
#' @param a.dist A character string specifying the accrual distribution. Must be one of
#'   "pieceunif" (piecewise uniform) or "piecepois" (piecewise Poisson). Default is "pieceunif".
#' @param a.time A numeric vector of time points defining the accrual intervals.
#'   For pieceunif: length should be one more than length of intensity/proportion.
#'   For piecepois: length should be equal to length of rate.
#'   Time points must be in increasing order.
#' @param intensity A numeric vector of accrual intensities (rate per time unit) for each interval
#'   when using piecewise uniform distribution. Cannot be used with proportion. Default is NULL.
#' @param proportion A numeric vector of accrual proportions for each interval that sum to 1
#'   when using piecewise uniform distribution. Cannot be used with intensity. Default is NULL.
#' @param rate A numeric vector of rates for each interval when using piecewise Poisson distribution.
#'   Length must equal length of a.time. All values must be non-negative. Default is NULL.
#' @param continuous A logical value for piecewise Poisson distribution. If TRUE, uses continuous
#'   Poisson process; if FALSE, uses discrete intervals. Default is FALSE.
#' @param e.time A named list where each element corresponds to a group, containing time points
#'   defining the survival time intervals. Last element must be Inf.
#' @param e.hazard A named list where each element corresponds to a group, containing hazard rates
#'   for survival time. All values must be positive.
#' @param d.time A named list where each element corresponds to a group, containing time points
#'   defining the dropout time intervals. Last element must be Inf. Can be NULL if no dropout
#'   is assumed. Default is NULL.
#' @param d.hazard A named list where each element corresponds to a group, containing hazard rates
#'   for dropout time. All values must be positive. Can be NULL if no dropout is assumed.
#'   Default is NULL.
#' @param seed A positive integer for random seed setting. Default is NULL.
#'
#' @return A tibble containing the simulated clinical trial data with the following columns:
#' \describe{
#'   \item{simID}{Simulation iteration ID (1 to nsim)}
#'   \item{patientID}{Patient ID within each simulation (1 to total sample size)}
#'   \item{group}{Group indicator (control, treatment, etc.)}
#'   \item{accrual.time}{Patient accrual time from study start}
#'   \item{survival.time}{Survival time from patient entry}
#'   \item{dropout.time}{Time to dropout from patient entry}
#'   \item{tte}{Time-to-event (minimum of survival time and dropout time)}
#'   \item{total}{Total time from study start (accrual time + time-to-event)}
#'   \item{dropout}{Dropout indicator (1 = dropout occurred, 0 = event occurred)}
#' }
#'
#' @details
#' This function uses group-wise operations to ensure statistical correctness. Each simulation
#' generates independent random numbers, avoiding issues with cumulative effects that can
#' occur with batch generation approaches.
#'
#' The function supports:
#' \describe{
#'   \item{Multiple Groups}{Any number of treatment groups can be specified.}
#'   \item{Flexible Accrual}{Choose between piecewise uniform or piecewise Poisson distributions.}
#'   \item{Optional Dropout}{Dropout can be disabled by setting d.time and d.hazard to NULL.}
#'   \item{Statistical Correctness}{Each simulation generates independent samples.}
#' }
#'
#' @examples
#' library(dplyr)
#'
#' # Example 1: Two-group trial with piecewise uniform accrual
#' trial1 <- simTrialData(
#'   nsim = 100,
#'   n = list(control = 250, treatment = 250),
#'   a.dist = 'pieceunif',
#'   a.time = c(0, 5, 10, 15, 20, 25),
#'   proportion = c(0.05, 0.1, 0.15, 0.3, 0.4),
#'   e.time = list(control = c(0, Inf), treatment = c(0, Inf)),
#'   e.hazard = list(control = log(2) / 4.3, treatment = log(2) / 5.811),
#'   d.time = list(control = c(0, Inf), treatment = c(0, Inf)),
#'   d.hazard = list(control = -log(1 - 0.1) / 12, treatment = -log(1 - 0.1) / 12),
#'   seed = 1
#' )
#'
#' # Example 2: Three-group trial without dropout
#' trial2 <- simTrialData(
#'   nsim = 50,
#'   n = list(control = 100, treatment1 = 100, treatment2 = 100),
#'   a.dist = 'pieceunif',
#'   a.time = c(0, 12, 24),
#'   intensity = c(12.5, 12.5),
#'   e.time = list(
#'     control = c(0, Inf),
#'     treatment1 = c(0, Inf),
#'     treatment2 = c(0, Inf)
#'   ),
#'   e.hazard = list(
#'     control = 0.08,
#'     treatment1 = 0.06,
#'     treatment2 = 0.05
#'   ),
#'   d.time = NULL,  # No dropout
#'   d.hazard = NULL
#' )
#'
#' # Example 3: Single-group trial with Poisson accrual
#' trial3 <- simTrialData(
#'   nsim = 100,
#'   n = list(treatment = 200),
#'   a.dist = 'piecepois',
#'   a.time = c(0, 6, 12, 18),
#'   rate = c(5, 10, 8, 2),
#'   continuous = TRUE,
#'   e.time = list(treatment = c(0, Inf)),
#'   e.hazard = list(treatment = 0.1),
#'   seed = 123
#' )
#'
#' @seealso
#' \code{\link{rpieceunif}} for piecewise uniform distribution,
#' \code{\link{rpieceexp}} for piecewise exponential distribution,
#' \code{\link{rpiecepois}} for piecewise Poisson distribution
#'
#' @import dplyr
#' @importFrom tibble tibble
#' @export
simTrialData <- function(nsim = 1000, n, a.dist = 'pieceunif', a.time,
                         intensity = NULL, proportion = NULL, rate = NULL, continuous = FALSE,
                         e.time, e.hazard, d.time = NULL, d.hazard = NULL, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Input validation
  if (!is.list(n) || length(n) < 1) {
    stop("n must be a list with at least 1 group")
  }

  if (!a.dist %in% c('pieceunif', 'piecepois')) {
    stop("a.dist must be either 'pieceunif' or 'piecepois'")
  }

  # Validate accrual parameters
  if (a.dist == 'pieceunif') {
    if (is.null(intensity) && is.null(proportion)) {
      stop("Either intensity or proportion must be specified for pieceunif distribution")
    }
    if (!is.null(intensity) && !is.null(proportion)) {
      stop("Cannot specify both intensity and proportion for pieceunif distribution")
    }
    # Check length consistency
    if (!is.null(intensity) && length(a.time) != length(intensity) + 1) {
      stop("For pieceunif with intensity: length of a.time must be one more than length of intensity")
    }
    if (!is.null(proportion) && length(a.time) != length(proportion) + 1) {
      stop("For pieceunif with proportion: length of a.time must be one more than length of proportion")
    }
  } else if (a.dist == 'piecepois') {
    if (is.null(rate)) {
      stop("rate must be specified for piecepois distribution")
    }
    if (length(a.time) != length(rate)) {
      stop("For piecepois distribution: length of a.time must equal length of rate")
    }
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

  # Calculate total sample size
  N <- sum(unlist(n))

  # Generate base dataset structure
  dataset <- tibble(
    simID = rep(1:nsim, each = N),
    patientID = rep(seq_len(N), nsim),
    group = rep(rep(names(n), unlist(n)), nsim)
  ) %>%
    group_by(simID) %>%
    mutate(
      # Generate accrual times
      accrual.time = if (a.dist == 'pieceunif') {
        rpieceunif(N, a.time, intensity, proportion)
      } else if (a.dist == 'piecepois') {
        rpiecepois(N, a.time, rate, continuous)
      },
      # Generate survival times
      survival.time = do.call(c, lapply(names(n), function(i) {
        rpieceexp(n[[i]], e.time[[i]], e.hazard[[i]])
      })),
      # Generate dropout times based on whether dropout is assumed
      dropout.time = if (has_dropout) {
        do.call(c, lapply(names(n), function(i) {
          rpieceexp(n[[i]], d.time[[i]], d.hazard[[i]])
        }))
      } else {
        rep(Inf, N)
      }
    ) %>%
    ungroup() %>%
    mutate(
      # Calculate time-to-event and total time
      tte = pmin(survival.time, dropout.time),
      total = accrual.time + tte,
      dropout = as.numeric(dropout.time < survival.time)
    )

  return(dataset)
}
