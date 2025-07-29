#' Simulate clinical trial data for multiple groups with optional subgroups using tidyverse
#'
#' This function generates simulation datasets for clinical trials with multiple treatment
#' groups and optional subgroups using tidyverse packages for efficient data manipulation.
#' It extends the functionality of simData to handle complex trial designs with time-to-event
#' endpoints, accounting for piecewise uniform distribution for patient accrual and piecewise
#' exponential distributions for patient survival time and dropout time.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 1000.
#' @param N A list where each element corresponds to a group, and if subgroups exist,
#'   each element is a named vector with subgroup sample sizes. If no subgroups,
#'   a vector of sample sizes for each group.
#' @param a.time A numeric vector of time points defining the accrual intervals for common accrual,
#'   OR a named list where each element corresponds to a subgroup with time points defining
#'   the accrual intervals. Time points must be in increasing order.
#' @param intensity A numeric vector of accrual intensities (rate per time unit) for each interval
#'   when using common accrual. Cannot be used with proportion or when a.time is a list. Default is NULL.
#' @param proportion A numeric vector of accrual proportions for each interval that sum to 1
#'   when using common accrual. Cannot be used with intensity or when a.time is a list. Default is NULL.
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
#' @param overall.time A numeric vector of overall study timeline when a.time is a named list with
#'   different subgroup timings. Required for delayed subgroup scenarios. Default is NULL.
#' @param overall.intensity A numeric vector of overall intensity targets when a.time is a named list.
#'   Used for delayed subgroup scenarios to maintain overall recruitment targets. Default is NULL.
#' @param seed A positive integer for random seed setting. Default is NULL.
#'
#' @return A tibble containing the simulated clinical trial data with the following columns:
#' \describe{
#'   \item{simID}{Simulation iteration ID (1 to nsim)}
#'   \item{group}{Group indicator (1, 2, 3, ...)}
#'   \item{subgroup}{Subgroup indicator (A, B, C, ...) - only present if subgroups are specified}
#'   \item{accrual.time}{Patient accrual time from study start}
#'   \item{surv.time}{Survival time from patient entry}
#'   \item{dropout.time}{Time to dropout from patient entry}
#'   \item{tte}{Time-to-event (minimum of survival time and dropout time)}
#'   \item{total}{Total time from study start (accrual time + time-to-event)}
#'   \item{dropout}{Dropout indicator (1 = dropout occurred, 0 = event occurred)}
#' }
#'
#' @details
#' This function is the main workhorse for clinical trial simulation in the FastSurvival
#' package. It supports:
#'
#' \describe{
#'   \item{Multiple Treatment Groups}{
#'     Each group can have different survival distributions and dropout rates.
#'   }
#'   \item{Subgroup Analysis}{
#'     Within each treatment group, subgroups can be defined with distinct
#'     characteristics (e.g., different biomarker status).
#'   }
#'   \item{Flexible Accrual Patterns}{
#'     Common accrual patterns across all groups/subgroups, OR subgroup-specific
#'     accrual patterns with delayed starts.
#'   }
#'   \item{Complex Survival Scenarios}{
#'     Supports delayed treatment effects, time-varying hazards, and different
#'     dropout patterns across groups. Dropout can be disabled by setting
#'     d.time and d.hazard to NULL.
#'   }
#' }
#'
#' The function uses piecewise exponential and uniform distributions to provide
#' maximum flexibility while maintaining computational efficiency through optimized
#' tidyverse operations.
#'
#' @examples
#' library(dplyr)
#'
#' # Example 1: Standard two-group trial without subgroups
#' trial1 <- simTrial(
#'   nsim = 100,
#'   N = list(control = 100, treatment = 100),
#'   a.time = c(0, 24),
#'   intensity = 200/24,
#'   e.time = list(control = c(0, Inf), treatment = c(0, Inf)),
#'   e.hazard = list(control = 0.08, treatment = 0.05),
#'   d.time = list(control = c(0, Inf), treatment = c(0, Inf)),
#'   d.hazard = list(control = 0.01, treatment = 0.01)
#' )
#'
#' # Example 2: Trial with subgroups and common accrual
#' trial2 <- simTrial(
#'   nsim = 100,
#'   N = list(
#'     control = c(A = 25, B = 112, C = 113),
#'     treatment = c(A = 25, B = 112, C = 113)
#'   ),
#'   a.time = c(0, 5, 10, 15, 20, 25),
#'   intensity = c(3.5, 14.3, 28.9, 43.6, 45),
#'   e.time = list(
#'     control = list(A = c(0, Inf), B = c(0, Inf), C = c(0, Inf)),
#'     treatment = list(A = c(0, Inf), B = c(0, Inf), C = c(0, Inf))
#'   ),
#'   e.hazard = list(
#'     control = list(A = log(2) / 5.811, B = log(2) / 5.811, C = log(2) / 5.811),
#'     treatment = list(A = log(2) / 4.3, B = log(2) / 4.3, C = log(2) / 4.3)
#'   ),
#'   d.time = list(
#'     control = list(A = c(0, Inf), B = c(0, Inf), C = c(0, Inf)),
#'     treatment = list(A = c(0, Inf), B = c(0, Inf), C = c(0, Inf))
#'   ),
#'   d.hazard = list(
#'     control = list(A = 0.01, B = 0.01, C = 0.01),
#'     treatment = list(A = 0.01, B = 0.01, C = 0.01)
#'   ),
#'   seed = 1
#' )
#'
#' # Example 3: Trial without dropout
#' trial3 <- simTrial(
#'   nsim = 100,
#'   N = list(control = 100, treatment = 100),
#'   a.time = c(0, 18),
#'   intensity = 200/18,
#'   e.time = list(control = c(0, Inf), treatment = c(0, Inf)),
#'   e.hazard = list(control = 0.08, treatment = 0.05),
#'   d.time = NULL,  # No dropout
#'   d.hazard = NULL
#' )
#'
#' # Example 4: Trial with delayed subgroup accrual
#' trial4 <- simTrial(
#'   nsim = 100,
#'   N = list(control = c(A = 50, B = 75), treatment = c(A = 50, B = 75)),
#'   a.time = list(A = c(6, 12, 24), B = c(0, 6, 12, 24)),  # Delayed start for A
#'   overall.time = c(0, 6, 12, 24),
#'   overall.intensity = c(10, 15, 20),
#'   e.time = list(
#'     control = list(A = c(0, Inf), B = c(0, Inf)),
#'     treatment = list(A = c(0, Inf), B = c(0, Inf))
#'   ),
#'   e.hazard = list(
#'     control = list(A = 0.08, B = 0.08),
#'     treatment = list(A = 0.05, B = 0.05)
#'   ),
#'   d.time = NULL,  # No dropout
#'   d.hazard = NULL
#' )
#'
#' @references
#' Luo, X., Mao, X., Chen, X., Qiu, J., Bai, S., & Quan, H. (2019).
#' Design and monitoring of survival trials in complex scenarios.
#' Statistics in Medicine, 38(2), 192-209.
#'
#' Friedman, L. M., Furberg, C. D., & DeMets, D. L. (2010).
#' Fundamentals of Clinical Trials (4th ed.). Springer.
#'
#' Proschan, M. A., Lan, K. K. G., & Wittes, J. T. (2006).
#' Statistical Monitoring of Clinical Trials: A Unified Approach. Springer.
#'
#' @seealso
#' \code{\link{simData}} for single-group simulation,
#' \code{\link{analysisData}} for creating analysis datasets,
#' \code{\link{lrtest}} and \code{\link{esthr}} for statistical analysis
#'
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom stats runif
#' @export
simTrial <- function(nsim = 1000, N, a.time, intensity = NULL, proportion = NULL,
                     e.time, e.hazard, d.time = NULL, d.hazard = NULL,
                     overall.time = NULL, overall.intensity = NULL, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Input validation
  if (!is.list(N) || length(N) < 1) {
    stop("N must be a list with at least 1 group")
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
  has_subgroups <- any(sapply(N, function(x) is.vector(x) && !is.null(names(x)) && length(x) > 1))

  # Check if we have delayed accrual (a.time is a list)
  has_delayed_accrual <- is.list(a.time)

  # Validate accrual parameters
  if (has_delayed_accrual) {
    # Delayed accrual scenario
    if (is.null(overall.time) || is.null(overall.intensity)) {
      stop("overall.time and overall.intensity must be provided when a.time is a list")
    }
    if (!is.null(intensity) || !is.null(proportion)) {
      stop("Cannot specify intensity or proportion when a.time is a list")
    }
  } else {
    # Common accrual scenario
    if (is.null(intensity) && is.null(proportion)) {
      stop("Either intensity or proportion must be specified for common accrual")
    }
    if (!is.null(intensity) && !is.null(proportion)) {
      stop("Cannot specify both intensity and proportion arguments")
    }
  }

  # Create parameter combinations
  param_combinations <- create_parameter_combinations(N, has_subgroups)

  # Generate data for all combinations
  all_data <- list()

  if (has_delayed_accrual) {
    # Use extractIntensity to calculate proportions for delayed accrual
    proportions_result <- extractIntensity(N, overall.time, overall.intensity, a.time)

    for (i in seq_len(nrow(param_combinations))) {
      combo <- param_combinations[i, ]

      # Extract sample size for this combination
      n_val <- extract_sample_size(N, combo$group, combo$subgroup, has_subgroups)

      if (has_subgroups && !is.na(combo$subgroup)) {
        # Get proportion for this subgroup from extractIntensity
        prop_col_name <- paste0(combo$subgroup, "_proportion")
        if (prop_col_name %in% names(proportions_result)) {
          subgroup_proportions <- proportions_result[[prop_col_name]]
        } else {
          subgroup_proportions <- rep(0, nrow(proportions_result))
        }

        # Calculate actual patient numbers for each period (using your approach)
        n_each_period <- pmax(1, round(subgroup_proportions * n_val))

        # Adjust the last period to ensure total equals n_val
        if (sum(n_each_period) != n_val) {
          last_active_idx <- max(which(n_each_period > 0))
          n_each_period[last_active_idx] <- n_each_period[last_active_idx] + (n_val - sum(n_each_period))
        }

        # Convert back to proportions, keeping only active periods
        active_periods <- n_each_period > 0
        final_proportions <- n_each_period[active_periods] / n_val

        # Get corresponding time points for active periods
        time_intervals <- sort(unique(c(overall.time, unlist(a.time))))
        active_time_start <- time_intervals[which(active_periods)]
        active_time_end <- time_intervals[which(active_periods) + 1]
        active_time_points <- c(active_time_start[1], active_time_end)

        # Generate accrual times using the calculated proportions
        accrual_times <- rpieceunif(
          n = n_val * nsim,
          time = active_time_points,
          proportion = final_proportions
        )
      } else {
        # Fallback for cases without subgroups
        accrual_times <- rpieceunif(
          n = n_val * nsim,
          time = c(min(overall.time), max(overall.time)),
          proportion = 1
        )
      }

      # Extract survival parameters
      e_time_val <- extract_parameter(e.time, combo$group, combo$subgroup, has_subgroups)
      e_hazard_val <- extract_parameter(e.hazard, combo$group, combo$subgroup, has_subgroups)

      # Extract dropout parameters if applicable
      if (has_dropout) {
        d_time_val <- extract_parameter(d.time, combo$group, combo$subgroup, has_subgroups)
        d_hazard_val <- extract_parameter(d.hazard, combo$group, combo$subgroup, has_subgroups)
        dropout_times <- rpieceexp(n_val * nsim, d_time_val, d_hazard_val)
      } else {
        dropout_times <- rep(Inf, n_val * nsim)
      }

      # Generate survival times
      surv_times <- rpieceexp(n_val * nsim, e_time_val, e_hazard_val)

      # Create group data using tibble
      group_data <- tibble(
        simID = rep(1:nsim, each = n_val),
        group = combo$group,
        accrual.time = accrual_times,
        surv.time = surv_times,
        dropout.time = dropout_times
      )

      # Add subgroup information if applicable
      if (has_subgroups) {
        group_data <- group_data %>%
          mutate(subgroup = combo$subgroup)
      }

      all_data[[i]] <- group_data
    }
  } else {
    # Common accrual scenario (existing logic)
    for (i in seq_len(nrow(param_combinations))) {
      combo <- param_combinations[i, ]

      # Extract sample size for this combination
      n_val <- extract_sample_size(N, combo$group, combo$subgroup, has_subgroups)

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

      # Use common accrual
      accrual_times <- rpieceunif(
        n = n_val * nsim,
        time = a.time,
        intensity = intensity,
        proportion = proportion
      )

      # Generate survival and dropout times
      surv_times <- rpieceexp(n_val * nsim, e_time_val, e_hazard_val)

      if (has_dropout) {
        dropout_times <- rpieceexp(n_val * nsim, d_time_val, d_hazard_val)
      } else {
        dropout_times <- rep(Inf, n_val * nsim)
      }

      # Create group data using tibble
      group_data <- tibble(
        simID = rep(1:nsim, each = n_val),
        group = combo$group,
        accrual.time = accrual_times,
        surv.time = surv_times,
        dropout.time = dropout_times
      )

      # Add subgroup information if applicable
      if (has_subgroups) {
        group_data <- group_data %>%
          mutate(subgroup = combo$subgroup)
      }

      all_data[[i]] <- group_data
    }
  }

  # Combine all data using bind_rows
  final_data <- bind_rows(all_data) %>%
    mutate(
      tte = pmin(surv.time, dropout.time),
      dropout = as.numeric(dropout.time < surv.time),
      total = accrual.time + tte
    )

  # Set column order
  if (has_subgroups) {
    final_data <- final_data %>%
      select(simID, group, subgroup, accrual.time, surv.time, dropout.time, tte, total, dropout)
  } else {
    final_data <- final_data %>%
      select(simID, group, accrual.time, surv.time, dropout.time, tte, total, dropout)
  }

  return(final_data)
}

# Helper function: Create parameter combinations using dplyr
create_parameter_combinations <- function(N, has_subgroups) {
  combinations <- data.frame(
    group = integer(),
    group_name = character(),
    subgroup = character(),
    stringsAsFactors = FALSE
  )

  if (has_subgroups) {
    for (group_idx in seq_along(N)) {
      group_data <- N[[group_idx]]
      group_name <- names(N)[group_idx]

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
    for (group_idx in seq_along(N)) {
      group_name <- names(N)[group_idx]
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
extract_sample_size <- function(N, group, subgroup, has_subgroups) {
  group_data <- N[[group]]

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
