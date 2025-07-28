#' Simulate clinical trial data for multiple groups with optional subgroups
#'
#' This function generates simulation datasets for clinical trials with multiple treatment
#' groups and optional subgroups using tidyverse packages. It provides a simplified
#' approach for complex trial designs with time-to-event endpoints.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 1000.
#' @param n A list where each element corresponds to a group, and if subgroups exist,
#'   each element is a named vector with subgroup sample sizes. If no subgroups,
#'   a numeric value for each group.
#' @param a.time A numeric vector of time points defining the accrual intervals.
#'   Time points must be in increasing order.
#' @param intensity A numeric vector of accrual intensities (rate per time unit) for each interval.
#'   Cannot be used with proportion. All values must be positive. Default is NULL.
#' @param proportion A numeric vector of accrual proportions for each interval that sum to 1.
#'   Cannot be used with intensity. All values must be non-negative and sum to 1. Default is NULL.
#' @param e.time A list where each element corresponds to a group (and subgroup if applicable),
#'   containing time points defining the survival time intervals. Last element must be Inf.
#' @param e.hazard A list where each element corresponds to a group (and subgroup if applicable),
#'   containing hazard rates for survival time. All values must be positive.
#' @param d.time A list where each element corresponds to a group (and subgroup if applicable),
#'   containing time points defining the dropout time intervals. Last element must be Inf.
#'   Can be NULL if no dropout is assumed. Default is NULL.
#' @param d.hazard A list where each element corresponds to a group (and subgroup if applicable),
#'   containing hazard rates for dropout time. All values must be positive.
#'   Can be NULL if no dropout is assumed. Default is NULL.
#' @param seed A positive integer for random seed setting. Default is NULL.
#'
#' @return A tibble containing the simulated clinical trial data with the following columns:
#' \describe{
#'   \item{simID}{Simulation iteration ID (1 to nsim)}
#'   \item{group}{Group indicator (1, 2, 3, ...)}
#'   \item{subgroup}{Subgroup indicator (A, B, C, ...) - only present if subgroups are specified}
#'   \item{patientID}{Patient ID within each group/subgroup combination}
#'   \item{accrual.time}{Patient accrual time from study start}
#'   \item{surv.time}{Survival time from patient entry}
#'   \item{dropout.time}{Time to dropout from patient entry}
#'   \item{tte}{Time-to-event (minimum of survival time and dropout time)}
#'   \item{total}{Total time from study start (accrual time + time-to-event)}
#'   \item{dropout}{Dropout indicator (1 = dropout occurred, 0 = event occurred)}
#' }
#'
#' @examples
#' library(dplyr)
#'
#' # Example 1: Standard two-group trial without subgroups
#' trial1 <- simTrial(
#'   nsim = 100,
#'   n = list(control = 100, treatment = 100),
#'   a.time = c(0, 24),
#'   intensity = 200/24,
#'   e.time = list(control = c(0, Inf), treatment = c(0, Inf)),
#'   e.hazard = list(control = 0.08, treatment = 0.05),
#'   d.time = list(control = c(0, Inf), treatment = c(0, Inf)),
#'   d.hazard = list(control = 0.01, treatment = 0.01)
#' )
#'
#' # Example 2: Trial with subgroups
#' trial2 <- simTrial(
#'   nsim = 100,
#'   n = list(
#'     control = c(A = 25, B = 75),
#'     treatment = c(A = 25, B = 75)
#'   ),
#'   a.time = c(0, 12, 24),
#'   intensity = c(10, 15),
#'   e.time = list(
#'     control = list(A = c(0, Inf), B = c(0, Inf)),
#'     treatment = list(A = c(0, Inf), B = c(0, Inf))
#'   ),
#'   e.hazard = list(
#'     control = list(A = 0.08, B = 0.08),
#'     treatment = list(A = 0.05, B = 0.05)
#'   ),
#'   d.time = NULL,
#'   d.hazard = NULL
#' )
#'
#' @seealso
#' \code{\link{simData}} for single-group simulation,
#' \code{\link{rpieceexp}} and \code{\link{rpieceunif}} for underlying distributions
#'
#' @import dplyr
#' @export
simTrial <- function(nsim = 1000, n, a.time, intensity = NULL, proportion = NULL,
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

  # Check if we have subgroups
  has_subgroups <- any(sapply(n, function(x) is.vector(x) && !is.null(names(x)) && length(x) > 1))

  # Create parameter combinations
  param_combinations <- create_parameter_combinations(n, has_subgroups)

  # Generate data for all combinations
  all_data <- list()

  for (i in seq_len(nrow(param_combinations))) {
    combo <- param_combinations[i, ]

    # Extract sample size for this combination
    n_val <- extract_sample_size(n, combo$group, combo$subgroup, has_subgroups)

    # Extract survival parameters
    e_time_val <- extract_parameter(e.time, combo$group, combo$subgroup, has_subgroups)
    e_hazard_val <- extract_parameter(e.hazard, combo$group, combo$subgroup, has_subgroups)

    # Extract dropout parameters if applicable
    if (has_dropout) {
      d_time_val <- extract_parameter(d.time, combo$group, combo$subgroup, has_subgroups)
      d_hazard_val <- extract_parameter(d.hazard, combo$group, combo$subgroup, has_subgroups)
    } else {
      d_time_val <- NULL
      d_hazard_val <- NULL
    }

    # Generate data using simData
    group_data <- simData(
      nsim = nsim,
      n = n_val,
      a.time = a.time,
      intensity = intensity,
      proportion = proportion,
      e.time = e_time_val,
      e.hazard = e_hazard_val,
      d.time = d_time_val,
      d.hazard = d_hazard_val
    )

    # Add group and subgroup information
    group_data <- group_data %>%
      mutate(
        group = combo$group,
        .before = patientID
      )

    if (has_subgroups) {
      group_data <- group_data %>%
        mutate(
          subgroup = combo$subgroup,
          .after = group
        )
    }

    all_data[[i]] <- group_data
  }

  # Combine all data
  final_data <- bind_rows(all_data)

  return(final_data)
}

# Helper function: Create parameter combinations
create_parameter_combinations <- function(n, has_subgroups) {
  combinations <- data.frame(
    group = integer(),
    group_name = character(),
    subgroup = character(),
    stringsAsFactors = FALSE
  )

  if (has_subgroups) {
    for (group_idx in seq_along(n)) {
      group_data <- n[[group_idx]]
      group_name <- names(n)[group_idx]

      if (is.vector(group_data) && !is.null(names(group_data)) && !is.list(group_data)) {
        # Named vector with subgroups
        subgroup_names <- names(group_data)
        for (subgroup_name in subgroup_names) {
          combinations <- rbind(combinations, data.frame(
            group = group_idx,
            group_name = group_name,
            subgroup = subgroup_name,
            stringsAsFactors = FALSE
          ))
        }
      } else {
        # Single value - treat as no subgroups for this group
        combinations <- rbind(combinations, data.frame(
          group = group_idx,
          group_name = group_name,
          subgroup = NA_character_,
          stringsAsFactors = FALSE
        ))
      }
    }
  } else {
    for (group_idx in seq_along(n)) {
      group_name <- names(n)[group_idx]
      combinations <- rbind(combinations, data.frame(
        group = group_idx,
        group_name = group_name,
        subgroup = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
  }

  return(combinations)
}

# Helper function: Extract sample size for specific combinations
extract_sample_size <- function(n, group, subgroup, has_subgroups) {
  group_data <- n[[group]]

  if (has_subgroups && !is.na(subgroup)) {
    if (is.vector(group_data) && !is.null(names(group_data)) && !is.list(group_data)) {
      # Named vector - extract by name
      if (subgroup %in% names(group_data)) {
        return(as.numeric(group_data[subgroup]))
      }
    }
  }

  # Return group level value
  if (is.vector(group_data) && length(group_data) == 1 && !is.list(group_data)) {
    return(as.numeric(group_data))
  } else if (is.list(group_data) && length(group_data) == 1) {
    return(as.numeric(group_data[[1]]))
  }

  return(as.numeric(group_data))
}

# Helper function: Extract parameter values for specific combinations
extract_parameter <- function(param_list, group, subgroup, has_subgroups) {
  group_data <- param_list[[group]]

  if (has_subgroups && !is.na(subgroup)) {
    if (is.list(group_data) && subgroup %in% names(group_data)) {
      return(group_data[[subgroup]])
    }
  }

  # Return group level parameter
  if (is.list(group_data) && length(group_data) == 1) {
    return(group_data[[1]])
  }

  return(group_data)
}
