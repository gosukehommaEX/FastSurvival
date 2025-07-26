#' Simulate clinical trial data for multiple groups with optional subgroups using data.table
#'
#' This function generates simulation datasets for clinical trials with multiple treatment
#' groups and optional subgroups using data.table for maximum performance and memory efficiency.
#' It extends the functionality of simData to handle complex trial designs with time-to-event
#' endpoints, accounting for piecewise uniform distribution for patient accrual and piecewise
#' exponential distributions for patient survival time and dropout time. The function supports
#' different accrual periods for different subgroups within treatment groups.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 1000.
#' @param N A list where each element corresponds to a group, and if subgroups exist,
#'   each element is a named vector with subgroup sample sizes. If no subgroups,
#'   a vector of sample sizes for each group.
#' @param a.time A numeric vector or named list of time points defining the accrual intervals.
#'   If numeric vector: common accrual periods for all groups/subgroups.
#'   If named list: subgroup-specific accrual periods allowing for delayed starts
#'   (e.g., list(A = c(6, 12, 24), B = c(0, 6, 12, 24))).
#'   Time points must be in increasing order.
#' @param intensity A numeric vector of accrual intensities (rate per time unit) for each interval.
#'   Cannot be used with proportion. Applied to all groups when a.time is numeric.
#'   Default is NULL.
#' @param proportion A numeric vector of accrual proportions for each interval that sum to 1.
#'   Cannot be used with intensity. Applied to all groups when a.time is numeric.
#'   Default is NULL.
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
#' @return A data.table containing the simulated clinical trial data with the following columns:
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
#'     Common accrual patterns across all groups, or different accrual periods
#'     for different subgroups to model delayed enrollment scenarios. Subgroup-specific
#'     timing applies to all treatment groups equally (randomization occurs after enrollment).
#'   }
#'   \item{Complex Survival Scenarios}{
#'     Supports delayed treatment effects, time-varying hazards, and different
#'     dropout patterns across groups. Dropout can be disabled by setting
#'     d.time and d.hazard to NULL.
#'   }
#' }
#'
#' The function uses piecewise exponential and uniform distributions to provide
#' maximum flexibility while maintaining computational efficiency.
#'
#' @examples
#' library(data.table)
#'
#' # Example 1: Standard two-group trial (existing functionality)
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
#' # Example 2: Trial with delayed subgroup accrual
#' trial2 <- simTrial(
#'   nsim = 100,
#'   N = list(control = c(A = 50, B = 75), treatment = c(A = 50, B = 75)),
#'   a.time = list(A = c(6, 12, 24), B = c(0, 6, 12, 24)),
#'   overall.time = c(0, 6, 12, 24),
#'   overall.intensity = c(10, 15, 20),
#'   e.time = list(
#'     control = list(A = c(0, Inf), B = c(0, Inf)),
#'     treatment = list(A = c(0, Inf), B = c(0, Inf))
#'   ),
#'   e.hazard = list(
#'     control = list(A = 0.08, B = 0.08),
#'     treatment = list(A = 0.05, B = 0.05)
#'   )
#' )
#'
#' @seealso
#' \code{\link{simData}} for single-group simulation,
#' \code{\link{analysisData}} for creating analysis datasets,
#' \code{\link{lrtest}} and \code{\link{esthr}} for statistical analysis
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
#' @import data.table
#' @export
simTrial <- function(nsim = 1e+3, N, a.time, intensity = NULL, proportion = NULL,
                     e.time, e.hazard, d.time = NULL, d.hazard = NULL,
                     overall.time = NULL, overall.intensity = NULL, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate input parameters
  if (is.null(intensity) && is.null(proportion) && is.null(overall.intensity)) {
    stop("Either intensity, proportion, or overall.intensity must be specified")
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

  # Detect delayed subgroup scenario
  delayed_scenario <- detect_delayed_subgroups(a.time, overall.time, overall.intensity)

  if (delayed_scenario$has_delays) {
    # Handle delayed subgroup accrual
    param_dt <- create_delayed_parameter_table(N, a.time, overall.time, overall.intensity,
                                               e.time, e.hazard, d.time, d.hazard, has_dropout)
  } else {
    # Use existing logic for uniform timing
    param_dt <- create_parameter_table_modified(N, a.time, intensity, proportion,
                                                e.time, e.hazard, d.time, d.hazard,
                                                delayed_scenario$has_subgroups, has_dropout)
  }

  # Generate data for all combinations using data.table operations
  all_data_list <- vector("list", nrow(param_dt))

  for (i in seq_len(nrow(param_dt))) {
    row <- param_dt[i]

    # Generate data using simData for this combination
    group_data <- simData(
      nsim = nsim,
      N = row$n_val[[1]],
      a.time = row$a_time_val[[1]],
      intensity = row$intensity_val[[1]],
      proportion = row$proportion_val[[1]],
      e.time = row$e_time_val[[1]],
      e.hazard = row$e_hazard_val[[1]],
      d.time = row$d_time_val[[1]],
      d.hazard = row$d_hazard_val[[1]]
    )

    # Add group and subgroup information using data.table operations
    if (!is.data.table(group_data)) {
      group_data <- as.data.table(group_data)
    }

    group_data[, group := row$group]
    if (delayed_scenario$has_subgroups) {
      group_data[, subgroup := row$subgroup]
    }

    all_data_list[[i]] <- group_data
  }

  # Combine all data using data.table's efficient rbindlist
  final_data <- rbindlist(all_data_list, use.names = TRUE, fill = TRUE)

  # Set optimal column order using data.table
  if (delayed_scenario$has_subgroups) {
    column_order <- c("simID", "group", "subgroup", "accrual.time",
                      "surv.time", "dropout.time", "tte", "total", "dropout")
  } else {
    column_order <- c("simID", "group", "accrual.time", "surv.time",
                      "dropout.time", "tte", "total", "dropout")
  }

  setcolorder(final_data, column_order)

  return(final_data)
}

# Helper function: Detect delayed subgroup scenario
detect_delayed_subgroups <- function(a.time, overall.time, overall.intensity) {

  # Check if delayed subgroup scenario
  has_delays <- FALSE
  has_subgroups <- FALSE

  if (is.list(a.time) && !is.null(overall.time) && !is.null(overall.intensity)) {
    # Check if a.time is a named list with subgroup timing
    if (!is.null(names(a.time))) {
      has_subgroups <- TRUE

      # Check for different start times across subgroups
      all_start_times <- sapply(a.time, function(x) x[1])
      has_delays <- length(unique(all_start_times)) > 1
    }
  }

  return(list(has_delays = has_delays, has_subgroups = has_subgroups))
}

# Helper function: Create parameter table for delayed scenarios
create_delayed_parameter_table <- function(N, a.time, overall.time, overall.intensity,
                                           e.time, e.hazard, d.time, d.hazard, has_dropout) {

  # Calculate proportion-based allocation
  allocation_result <- calculate_delayed_allocation(N, a.time, overall.time, overall.intensity)

  # Extract combinations
  combinations <- extract_combinations_dt(N, allocation_result$has_subgroups)

  # Create parameter data.table
  param_dt <- data.table(
    group = combinations$group,
    subgroup = combinations$subgroup,
    combination_key = combinations$combination_key
  )

  # Add parameter values using data.table operations
  param_dt[, `:=`(
    n_val = extract_parameter_values(N, combination_key, allocation_result$has_subgroups),
    a_time_val = extract_subgroup_a_time(a.time, combination_key, allocation_result$has_subgroups),
    e_time_val = extract_parameter_values(e.time, combination_key, allocation_result$has_subgroups),
    e_hazard_val = extract_parameter_values(e.hazard, combination_key, allocation_result$has_subgroups)
  )]

  # Handle dropout parameters
  if (has_dropout) {
    param_dt[, `:=`(
      d_time_val = extract_parameter_values(d.time, combination_key, allocation_result$has_subgroups),
      d_hazard_val = extract_parameter_values(d.hazard, combination_key, allocation_result$has_subgroups)
    )]
  } else {
    param_dt[, `:=`(
      d_time_val = list(NULL),
      d_hazard_val = list(NULL)
    )]
  }

  # Add calculated proportions
  param_dt[, `:=`(
    intensity_val = list(NULL),
    proportion_val = extract_subgroup_proportions(allocation_result$subgroup_proportions,
                                                  combination_key, allocation_result$has_subgroups)
  )]

  return(param_dt)
}

# Helper function: Calculate delayed allocation
calculate_delayed_allocation <- function(N, a.time, overall.time, overall.intensity) {

  # Convert overall intensity to proportions
  interval_lengths <- diff(overall.time)
  expected_counts <- overall.intensity * interval_lengths
  total_expected <- sum(expected_counts)
  overall_proportions <- expected_counts / total_expected

  # Calculate target allocation per interval
  total_N <- sum(unlist(N))
  target_allocation <- overall_proportions * total_N

  # Initialize subgroup allocation - now subgroup-specific only
  subgroup_allocation <- list()
  has_subgroups <- TRUE  # Always true for delayed scenario

  # Create allocation structure for each subgroup
  for (sg_name in names(a.time)) {
    sg_intervals <- length(a.time[[sg_name]]) - 1
    subgroup_allocation[[sg_name]] <- numeric(sg_intervals)
  }

  # Track remaining sample sizes by subgroup across all groups
  remaining_N_by_subgroup <- list()
  for (sg_name in names(a.time)) {
    total_sg_N <- sum(sapply(N, function(group) {
      if (sg_name %in% names(group)) group[[sg_name]] else 0
    }))
    remaining_N_by_subgroup[[sg_name]] <- total_sg_N
  }

  # Distribute allocation across intervals
  for (i in 1:length(target_allocation)) {
    target_for_interval <- target_allocation[i]

    # Find active subgroups and calculate remaining needs
    active_info <- get_active_subgroups_info_simplified(i, a.time, overall.time, remaining_N_by_subgroup)

    if (active_info$total_remaining <= 0) next

    # Allocate based on remaining needs
    if (active_info$total_remaining <= target_for_interval) {
      # Satisfy all remaining needs
      for (j in 1:length(active_info$active_keys)) {
        key_info <- active_info$active_keys[[j]]
        allocation_amount <- key_info$remaining

        subgroup_allocation[[key_info$subgroup]][key_info$interval_idx] <- allocation_amount
        remaining_N_by_subgroup[[key_info$subgroup]] <- 0
      }
    } else {
      # Proportional allocation
      for (j in 1:length(active_info$active_keys)) {
        key_info <- active_info$active_keys[[j]]
        proportion <- key_info$remaining / active_info$total_remaining
        allocation_amount <- target_for_interval * proportion

        subgroup_allocation[[key_info$subgroup]][key_info$interval_idx] <- allocation_amount
        remaining_N_by_subgroup[[key_info$subgroup]] <- remaining_N_by_subgroup[[key_info$subgroup]] - allocation_amount
      }
    }
  }

  # Convert to proportions
  subgroup_proportions <- convert_allocation_to_proportions_simplified(subgroup_allocation)

  return(list(
    subgroup_proportions = subgroup_proportions,
    has_subgroups = has_subgroups
  ))
}

# Helper function: Get active subgroups info for an interval (simplified for subgroup-only timing)
get_active_subgroups_info_simplified <- function(interval_idx, a.time, overall.time, remaining_N_by_subgroup) {

  interval_start <- overall.time[interval_idx]
  interval_end <- overall.time[interval_idx + 1]

  active_keys <- list()
  total_remaining <- 0

  for (sg_name in names(a.time)) {
    sg_times <- a.time[[sg_name]]
    sg_interval_idx <- find_overlapping_interval(interval_start, interval_end, sg_times)

    if (!is.na(sg_interval_idx)) {
      remaining <- remaining_N_by_subgroup[[sg_name]]
      if (remaining > 0) {
        active_keys <- append(active_keys, list(list(
          subgroup = sg_name,
          interval_idx = sg_interval_idx,
          remaining = remaining
        )), after = length(active_keys))
        total_remaining <- total_remaining + remaining
      }
    }
  }

  return(list(active_keys = active_keys, total_remaining = total_remaining))
}

# Helper function: Find overlapping interval index
find_overlapping_interval <- function(interval_start, interval_end, sg_times) {
  for (j in 1:(length(sg_times) - 1)) {
    sg_start <- sg_times[j]
    sg_end <- sg_times[j + 1]

    if (sg_start < interval_end && sg_end > interval_start) {
      return(j)
    }
  }
  return(NA)
}

# Helper function: Convert allocation to proportions (simplified for subgroup-only)
convert_allocation_to_proportions_simplified <- function(subgroup_allocation) {

  subgroup_proportions <- list()

  for (sg_name in names(subgroup_allocation)) {
    allocations <- subgroup_allocation[[sg_name]]
    total_allocated <- sum(allocations)
    if (total_allocated > 0) {
      subgroup_proportions[[sg_name]] <- allocations / total_allocated
    } else {
      subgroup_proportions[[sg_name]] <- allocations
    }
  }

  return(subgroup_proportions)
}

# Helper function: Extract subgroup-specific a.time
extract_subgroup_a_time <- function(a.time, combination_keys, has_subgroups) {

  result <- vector("list", length(combination_keys))

  for (i in seq_along(combination_keys)) {
    key <- combination_keys[i]

    if (has_subgroups) {
      # Parse subgroup from combination key
      parts <- strsplit(key, "_")[[1]]
      subgroup_name <- parts[2]

      # For delayed scenario, a.time is subgroup-specific
      result[[i]] <- a.time[[subgroup_name]]
    } else {
      # This shouldn't happen in delayed scenario, but keep for safety
      result[[i]] <- a.time
    }
  }

  return(result)
}
# Helper function: Extract subgroup proportions (simplified for subgroup-only)
extract_subgroup_proportions <- function(subgroup_proportions, combination_keys, has_subgroups) {

  result <- vector("list", length(combination_keys))

  for (i in seq_along(combination_keys)) {
    key <- combination_keys[i]

    if (has_subgroups) {
      # Parse subgroup from combination key
      parts <- strsplit(key, "_")[[1]]
      subgroup_name <- parts[2]

      result[[i]] <- subgroup_proportions[[subgroup_name]]
    } else {
      # This shouldn't happen in delayed scenario, but keep for safety
      group_idx <- as.integer(gsub("group", "", key))
      group_names <- names(subgroup_proportions)
      group_name <- group_names[group_idx]

      result[[i]] <- subgroup_proportions[[group_name]]
    }
  }

  return(result)
}

# Reuse existing helper functions from original simTrial with modifications
create_parameter_table_modified <- function(N, a.time, intensity, proportion,
                                            e.time, e.hazard, d.time, d.hazard,
                                            has_subgroups, has_dropout) {

  # Extract combinations into data.table format
  combinations <- extract_combinations_dt(N, has_subgroups)

  # Create parameter data.table
  param_dt <- data.table(
    group = combinations$group,
    subgroup = combinations$subgroup,
    combination_key = combinations$combination_key
  )

  # Add parameter values using data.table operations
  param_dt[, `:=`(
    n_val = extract_parameter_values(N, combination_key, has_subgroups),
    a_time_val = list(a.time),  # Common accrual time for all groups
    e_time_val = extract_parameter_values(e.time, combination_key, has_subgroups),
    e_hazard_val = extract_parameter_values(e.hazard, combination_key, has_subgroups)
  )]

  # Handle dropout parameters
  if (has_dropout) {
    param_dt[, `:=`(
      d_time_val = extract_parameter_values(d.time, combination_key, has_subgroups),
      d_hazard_val = extract_parameter_values(d.hazard, combination_key, has_subgroups)
    )]
  } else {
    param_dt[, `:=`(
      d_time_val = list(NULL),
      d_hazard_val = list(NULL)
    )]
  }

  # Handle accrual parameters (common for all groups)
  if (!is.null(intensity)) {
    param_dt[, `:=`(
      intensity_val = list(intensity),  # Common intensity for all groups
      proportion_val = list(NULL)
    )]
  } else {
    param_dt[, `:=`(
      intensity_val = list(NULL),
      proportion_val = list(proportion)  # Common proportion for all groups
    )]
  }

  return(param_dt)
}

# Internal function: Extract group and subgroup combinations using data.table
extract_combinations_dt <- function(param_list, has_subgroups) {

  group_vals <- integer(0)
  subgroup_vals <- character(0)
  combination_keys <- character(0)

  if (has_subgroups) {
    for (g in seq_along(param_list)) {
      if (is.list(param_list[[g]])) {
        subgroup_names <- names(param_list[[g]])
        if (is.null(subgroup_names)) {
          subgroup_names <- LETTERS[seq_along(param_list[[g]])]
        }
        for (s in seq_along(param_list[[g]])) {
          group_vals <- c(group_vals, g)
          subgroup_vals <- c(subgroup_vals, subgroup_names[s])
          combination_keys <- c(combination_keys, paste0("group", g, "_", subgroup_names[s]))
        }
      } else {
        subgroup_names <- names(param_list[[g]])
        if (is.null(subgroup_names)) {
          subgroup_names <- LETTERS[seq_along(param_list[[g]])]
        }
        for (s in seq_along(param_list[[g]])) {
          group_vals <- c(group_vals, g)
          subgroup_vals <- c(subgroup_vals, subgroup_names[s])
          combination_keys <- c(combination_keys, paste0("group", g, "_", subgroup_names[s]))
        }
      }
    }
  } else {
    for (g in seq_along(param_list)) {
      group_vals <- c(group_vals, g)
      subgroup_vals <- c(subgroup_vals, NA_character_)
      combination_keys <- c(combination_keys, paste0("group", g))
    }
  }

  return(list(
    group = group_vals,
    subgroup = subgroup_vals,
    combination_key = combination_keys
  ))
}

# Internal function: Extract parameter values for specific combinations
extract_parameter_values <- function(param_list, combination_keys, has_subgroups) {

  result <- vector("list", length(combination_keys))

  for (i in seq_along(combination_keys)) {
    key <- combination_keys[i]

    if (has_subgroups) {
      # Parse group and subgroup from key
      parts <- strsplit(key, "_")[[1]]
      group_idx <- as.integer(gsub("group", "", parts[1]))
      subgroup_name <- parts[2]

      if (is.list(param_list[[group_idx]])) {
        if (subgroup_name %in% names(param_list[[group_idx]])) {
          result[[i]] <- param_list[[group_idx]][[subgroup_name]]
        } else {
          subgroup_idx <- match(subgroup_name, LETTERS)
          result[[i]] <- param_list[[group_idx]][[subgroup_idx]]
        }
      } else {
        if (subgroup_name %in% names(param_list[[group_idx]])) {
          result[[i]] <- param_list[[group_idx]][subgroup_name]
        } else {
          subgroup_idx <- match(subgroup_name, LETTERS)
          result[[i]] <- param_list[[group_idx]][subgroup_idx]
        }
      }
    } else {
      group_idx <- as.integer(gsub("group", "", key))
      result[[i]] <- param_list[[group_idx]]
    }
  }

  return(result)
}
