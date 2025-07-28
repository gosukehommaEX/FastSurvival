#' Simulate clinical trial data for multiple groups with optional subgroups using tidyverse
#'
#' This function generates simulation datasets for clinical trials with multiple treatment
#' groups and optional subgroups using tidyverse packages for maximum performance and
#' flexibility. It extends the functionality of simData to handle complex trial designs
#' with time-to-event endpoints, accounting for piecewise uniform distribution for patient
#' accrual and piecewise exponential distributions for patient survival time and dropout time.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 1000.
#' @param n A list where each element corresponds to a group, and if subgroups exist,
#'   each element is a named vector with subgroup sample sizes. If no subgroups,
#'   a vector of sample sizes for each group.
#' @param a.time A numeric vector of time points defining the accrual intervals for common accrual,
#'   OR a named list where each element corresponds to a subgroup with time points defining
#'   the accrual intervals. Time points must be in increasing order.
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
#' @details
#' This function is the main workhorse for clinical trial simulation. It supports:
#'
#' \describe{
#'   \item{Multiple Treatment Groups}{
#'     Each group can have different survival distributions and dropout rates.
#'   }
#'   \item{Subgroup Analysis}{
#'     Within each treatment group, subgroups can be defined with distinct
#'     characteristics and different accrual timings.
#'   }
#'   \item{Flexible Accrual Patterns}{
#'     Common accrual patterns across all groups/subgroups, OR subgroup-specific
#'     accrual patterns with automatic intensity adjustment to maintain overall targets.
#'   }
#'   \item{Complex Survival Scenarios}{
#'     Supports delayed treatment effects, time-varying hazards, and different
#'     dropout patterns across groups. Dropout can be disabled by setting
#'     d.time and d.hazard to NULL.
#'   }
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
#' # Example 2: Trial with subgroups and common accrual
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
#' # Example 3: Trial with delayed subgroup accrual
#' trial3 <- simTrial(
#'   nsim = 50,
#'   n = list(control = c(A = 50, B = 50), treatment = c(A = 50, B = 50)),
#'   a.time = list(A = c(3, 5, 10, 15), B = c(0, 5, 10, 15)),
#'   intensity = c(10, 20, 30),
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
#' \code{\link{analysisData}} for creating analysis datasets,
#' \code{\link{lrtest}} and \code{\link{esthr}} for statistical analysis
#'
#' @import dplyr
#' @import tidyr
#' @importFrom tibble tibble
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
  has_subgroups <- check_subgroup_structure(n)

  # Create parameter combinations
  param_combinations <- create_parameter_combinations(n, has_subgroups)

  # Calculate accrual parameters for each combination
  accrual_params <- calculate_accrual_parameters(
    param_combinations, a.time, intensity, proportion, has_subgroups
  )

  # Generate data for all combinations
  all_data <- do.call(rbind, lapply(seq_len(nrow(param_combinations)), function(i) {
    combo <- param_combinations[i, ]
    accrual_param <- accrual_params[i, ]

    # Extract sample size for this combination
    n_val <- extract_parameter_value(n, combo$group, combo$subgroup, has_subgroups)

    # Extract survival parameters
    e_time_val <- extract_parameter_value(e.time, combo$group, combo$subgroup, has_subgroups)
    e_hazard_val <- extract_parameter_value(e.hazard, combo$group, combo$subgroup, has_subgroups)

    # Extract dropout parameters if applicable
    if (has_dropout) {
      d_time_val <- extract_parameter_value(d.time, combo$group, combo$subgroup, has_subgroups)
      d_hazard_val <- extract_parameter_value(d.hazard, combo$group, combo$subgroup, has_subgroups)
    } else {
      d_time_val <- NULL
      d_hazard_val <- NULL
    }

    # Extract accrual parameters properly - need double extraction for nested lists
    a_time_val <- accrual_param$a_time_adj[[1]]
    intensity_val <- accrual_param$intensity_adj[[1]]
    proportion_val <- accrual_param$proportion_adj[[1]]

    # If they are still lists (nested), extract again
    if (is.list(a_time_val)) a_time_val <- a_time_val[[1]]
    if (is.list(intensity_val)) intensity_val <- intensity_val[[1]]
    if (is.list(proportion_val)) proportion_val <- proportion_val[[1]]

    # Debug: Check the values and types after extraction (comment out for production)
    # cat("Group:", combo$group, "Subgroup:", combo$subgroup, "\n")
    # cat("a_time_val type:", class(a_time_val), "length:", length(a_time_val), "values:", paste(a_time_val, collapse = ", "), "\n")
    # cat("intensity_val type:", class(intensity_val), "length:", length(intensity_val), "values:", paste(intensity_val, collapse = ", "), "\n")
    # cat("proportion_val:", paste(proportion_val, collapse = ", "), "\n")

    # Ensure only one of intensity or proportion is non-NULL
    if (!is.null(intensity_val) && !is.null(proportion_val)) {
      # This shouldn't happen, but as a safety measure
      proportion_val <- NULL
    }

    # Validate parameters before calling simData
    if (!is.null(intensity_val)) {
      if (length(a_time_val) != length(intensity_val) + 1) {
        stop("Length mismatch: a.time length = ", length(a_time_val),
             ", intensity length = ", length(intensity_val))
      }
    }

    if (!is.null(proportion_val)) {
      if (length(a_time_val) != length(proportion_val) + 1) {
        stop("Length mismatch: a.time length = ", length(a_time_val),
             ", proportion length = ", length(proportion_val))
      }
    }

    # Generate data using simData
    group_data <- simData(
      nsim = nsim,
      n = n_val,
      a.time = a_time_val,
      intensity = intensity_val,
      proportion = proportion_val,
      e.time = e_time_val,
      e.hazard = e_hazard_val,
      d.time = d_time_val,
      d.hazard = d_hazard_val
    )

    # Add group and subgroup information
    result <- group_data %>%
      mutate(
        group = combo$group,
        .before = patientID
      )

    if (has_subgroups) {
      result <- result %>%
        mutate(
          subgroup = combo$subgroup,
          .after = group
        )
    }

    return(result)
  }))

  return(all_data)
}

# Helper function: Check if n has subgroup structure
check_subgroup_structure <- function(n) {
  if (!is.list(n)) return(FALSE)

  # Check if any element in n is a named vector (indicating subgroups)
  any(vapply(n, function(x) {
    is.vector(x) && !is.null(names(x)) && length(x) > 1 && !is.list(x)
  }, logical(1)))
}

# Helper function: Create parameter combinations
create_parameter_combinations <- function(n, has_subgroups) {
  combinations <- data.frame(stringsAsFactors = FALSE)

  if (has_subgroups) {
    for (group_idx in seq_along(n)) {
      group_data <- n[[group_idx]]
      group_name <- names(n)[group_idx]

      if (is.vector(group_data) && !is.null(names(group_data)) && !is.list(group_data)) {
        # Named vector with subgroups
        subgroup_names <- names(group_data)
        for (subgroup_name in subgroup_names) {
          new_row <- data.frame(
            group = group_idx,
            group_name = group_name,
            subgroup = subgroup_name,
            stringsAsFactors = FALSE
          )
          combinations <- rbind(combinations, new_row)
        }
      } else {
        # Single value - treat as no subgroups for this group
        new_row <- data.frame(
          group = group_idx,
          group_name = group_name,
          subgroup = NA_character_,
          stringsAsFactors = FALSE
        )
        combinations <- rbind(combinations, new_row)
      }
    }
  } else {
    for (group_idx in seq_along(n)) {
      group_name <- names(n)[group_idx]
      new_row <- data.frame(
        group = group_idx,
        group_name = group_name,
        subgroup = NA_character_,
        stringsAsFactors = FALSE
      )
      combinations <- rbind(combinations, new_row)
    }
  }

  return(combinations)
}

# Helper function: Calculate accrual parameters
calculate_accrual_parameters <- function(param_combinations, a.time, intensity, proportion, has_subgroups) {

  # Check if we have subgroup-specific accrual times
  subgroup_specific_accrual <- is.list(a.time) && has_subgroups

  if (subgroup_specific_accrual) {
    # Calculate overall time points and adjust intensities/proportions
    return(calculate_subgroup_accrual_adjustment(param_combinations, a.time, intensity, proportion))
  } else {
    # Common accrual for all combinations
    n_combinations <- nrow(param_combinations)

    # Ensure intensity/proportion vectors have correct length
    if (!is.null(intensity)) {
      intensity_vector <- if (length(intensity) == length(a.time) - 1) {
        intensity
      } else if (length(intensity) == 1) {
        rep(intensity, length(a.time) - 1)
      } else {
        stop("Intensity length must be 1 or equal to length(a.time) - 1")
      }
    } else {
      intensity_vector <- NULL
    }

    if (!is.null(proportion)) {
      proportion_vector <- if (length(proportion) == length(a.time) - 1) {
        proportion
      } else {
        stop("Proportion length must be equal to length(a.time) - 1")
      }
    } else {
      proportion_vector <- NULL
    }

    # Create result with proper list columns
    a_time_list <- replicate(n_combinations, list(a.time), simplify = FALSE)
    intensity_list <- replicate(n_combinations, list(intensity_vector), simplify = FALSE)
    proportion_list <- replicate(n_combinations, list(proportion_vector), simplify = FALSE)

    result <- data.frame(
      row_id = seq_len(n_combinations),
      stringsAsFactors = FALSE
    )
    result$a_time_adj <- a_time_list
    result$intensity_adj <- intensity_list
    result$proportion_adj <- proportion_list
    result$row_id <- NULL  # Remove temporary column

    return(result)
  }
}

# Helper function: Calculate subgroup accrual adjustment
calculate_subgroup_accrual_adjustment <- function(param_combinations, a.time, intensity, proportion) {

  # Calculate overall time points
  all_time_points <- unique(sort(unlist(a.time)))

  # Get unique subgroups
  subgroups <- unique(param_combinations$subgroup[!is.na(param_combinations$subgroup)])

  # Calculate interval information
  n_intervals <- length(all_time_points) - 1
  interval_lengths <- diff(all_time_points)

  # For each interval, determine which subgroups are active
  interval_subgroups <- lapply(seq_len(n_intervals), function(i) {
    interval_start <- all_time_points[i]
    interval_end <- all_time_points[i + 1]

    active_subgroups <- character(0)
    for (subgroup in subgroups) {
      if (subgroup %in% names(a.time)) {
        subgroup_times <- a.time[[subgroup]]
        subgroup_start <- min(subgroup_times)
        subgroup_end <- max(subgroup_times)

        # Check if this subgroup is active in this interval
        if (interval_start >= subgroup_start && interval_end <= subgroup_end) {
          active_subgroups <- c(active_subgroups, subgroup)
        }
      }
    }
    return(active_subgroups)
  })

  # Calculate adjusted intensities/proportions for each subgroup
  result <- do.call(rbind, lapply(seq_len(nrow(param_combinations)), function(i) {
    combo <- param_combinations[i, ]

    if (is.na(combo$subgroup)) {
      # No subgroup - use common accrual
      result <- data.frame(row_id = 1, stringsAsFactors = FALSE)
      result$a_time_adj <- list(all_time_points)
      result$intensity_adj <- list(if (!is.null(intensity)) intensity else NULL)
      result$proportion_adj <- list(if (!is.null(proportion)) proportion else NULL)
      result$row_id <- NULL
      return(result)
    }

    subgroup_name <- combo$subgroup

    if (!subgroup_name %in% names(a.time)) {
      # Subgroup not in a.time list - use common accrual
      result <- data.frame(row_id = 1, stringsAsFactors = FALSE)
      result$a_time_adj <- list(all_time_points)
      result$intensity_adj <- list(if (!is.null(intensity)) intensity else NULL)
      result$proportion_adj <- list(if (!is.null(proportion)) proportion else NULL)
      result$row_id <- NULL
      return(result)
    }

    # Get subgroup-specific time points
    subgroup_times <- a.time[[subgroup_name]]
    subgroup_start <- min(subgroup_times)
    subgroup_end <- max(subgroup_times)

    # Find intervals where this subgroup is active
    active_intervals <- which(
      all_time_points[-length(all_time_points)] >= subgroup_start &
        all_time_points[-1] <= subgroup_end
    )

    if (length(active_intervals) == 0) {
      # No active intervals - return zero intensity/proportion
      zero_intensity <- rep(0, n_intervals)
      zero_proportion <- rep(0, n_intervals)
      zero_proportion[1] <- 1  # Avoid sum != 1 error

      result <- data.frame(row_id = 1, stringsAsFactors = FALSE)
      result$a_time_adj <- list(all_time_points)
      result$intensity_adj <- list(if (!is.null(intensity)) zero_intensity else NULL)
      result$proportion_adj <- list(if (!is.null(proportion)) zero_proportion else NULL)
      result$row_id <- NULL
      return(result)
    }

    # Calculate adjusted parameters
    if (!is.null(intensity)) {
      # Adjust intensity based on number of active subgroups in each interval
      adj_intensity <- numeric(n_intervals)

      for (j in seq_len(n_intervals)) {
        if (j %in% active_intervals) {
          n_active_subgroups <- length(interval_subgroups[[j]])
          target_intensity <- intensity[min(j, length(intensity))]
          adj_intensity[j] <- target_intensity / max(1, n_active_subgroups)
        } else {
          adj_intensity[j] <- 0
        }
      }

      result <- data.frame(row_id = 1, stringsAsFactors = FALSE)
      result$a_time_adj <- list(all_time_points)
      result$intensity_adj <- list(adj_intensity)
      result$proportion_adj <- list(NULL)
      result$row_id <- NULL
      return(result)

    } else {
      # Adjust proportions
      adj_proportion <- numeric(n_intervals)
      total_target_proportion <- sum(proportion)

      active_proportion_sum <- sum(proportion[active_intervals])
      if (active_proportion_sum > 0) {
        for (j in active_intervals) {
          target_prop <- proportion[min(j, length(proportion))]
          n_active_subgroups <- length(interval_subgroups[[j]])
          adj_proportion[j] <- target_prop / max(1, n_active_subgroups)
        }
      }

      # Normalize to sum to 1
      if (sum(adj_proportion) > 0) {
        adj_proportion <- adj_proportion / sum(adj_proportion)
      } else {
        adj_proportion[1] <- 1  # Fallback
      }

      result <- data.frame(row_id = 1, stringsAsFactors = FALSE)
      result$a_time_adj <- list(all_time_points)
      result$intensity_adj <- list(NULL)
      result$proportion_adj <- list(adj_proportion)
      result$row_id <- NULL
      return(result)
    }
  }))

  return(result)
}

# Helper function: Extract parameter values for specific combinations
extract_parameter_value <- function(param_list, group, subgroup, has_subgroups) {

  if (has_subgroups && !is.na(subgroup)) {
    # Extract subgroup-specific parameter
    group_data <- param_list[[group]]

    if (is.vector(group_data) && !is.null(names(group_data)) && !is.list(group_data)) {
      # Named vector - extract by name (for n parameter)
      if (subgroup %in% names(group_data)) {
        return(as.numeric(unname(group_data[subgroup])))
      }
    } else if (is.list(group_data)) {
      # List structure (for e.time, e.hazard, d.time, d.hazard)
      if (subgroup %in% names(group_data)) {
        return(group_data[[subgroup]])
      }
    }

    # Fallback to group level if subgroup not found
    if (is.list(group_data) && length(group_data) == 1) {
      return(group_data[[1]])
    } else if (is.vector(group_data) && length(group_data) == 1) {
      return(as.numeric(group_data))
    }
  } else {
    # Extract group-level parameter
    group_data <- param_list[[group]]

    if (is.vector(group_data) && length(group_data) == 1 && !is.list(group_data)) {
      return(as.numeric(unname(group_data)))
    } else if (is.list(group_data) && length(group_data) == 1) {
      return(group_data[[1]])
    } else {
      return(group_data)
    }
  }

  return(NA)
}
