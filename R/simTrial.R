#' Simulate clinical trial data for multiple groups with optional subgroups using data.table
#'
#' This function generates simulation datasets for clinical trials with multiple treatment
#' groups and optional subgroups using data.table for maximum performance and memory efficiency.
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
#' maximum flexibility while maintaining computational efficiency.
#'
#' @examples
#' library(data.table)
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

  # Input validation
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

  # Check if we have subgroups in N parameter structure
  has_subgroups <- check_subgroup_structure(N)

  # Determine accrual scenario
  is_delayed_subgroup <- is.list(a.time) && !is.null(overall.time) && !is.null(overall.intensity)

  # Create parameter table based on scenario
  if (is_delayed_subgroup) {
    param_dt <- create_parameter_table_delayed(N, a.time, overall.time, overall.intensity,
                                               e.time, e.hazard, d.time, d.hazard,
                                               has_subgroups, has_dropout)
  } else {
    param_dt <- create_parameter_table_standard(N, a.time, intensity, proportion,
                                                e.time, e.hazard, d.time, d.hazard,
                                                has_subgroups, has_dropout)
  }

  # Generate data for all combinations using data.table operations
  all_data_list <- vector("list", nrow(param_dt))

  for (i in seq_len(nrow(param_dt))) {
    row <- param_dt[i]

    # Extract values
    n_val <- row$n_val[[1]]
    a_time_val <- row$a_time_val[[1]]
    intensity_val <- row$intensity_val[[1]]
    proportion_val <- row$proportion_val[[1]]
    e_time_val <- row$e_time_val[[1]]
    e_hazard_val <- row$e_hazard_val[[1]]
    d_time_val <- row$d_time_val[[1]]
    d_hazard_val <- row$d_hazard_val[[1]]

    # Ensure hazard values are numeric vectors, not lists
    if (is.list(e_hazard_val)) {
      e_hazard_val <- unlist(e_hazard_val)
    }
    if (!is.null(d_hazard_val) && is.list(d_hazard_val)) {
      d_hazard_val <- unlist(d_hazard_val)
    }

    # Generate data using simData for this combination
    group_data <- simData(
      nsim = nsim,
      N = n_val,
      a.time = a_time_val,
      intensity = intensity_val,
      proportion = proportion_val,
      e.time = e_time_val,
      e.hazard = e_hazard_val,
      d.time = d_time_val,
      d.hazard = d_hazard_val
    )

    # Add group and subgroup information using data.table operations
    if (!is.data.table(group_data)) {
      group_data <- as.data.table(group_data)
    }

    group_data[, group := row$group]
    if (has_subgroups) {
      group_data[, subgroup := row$subgroup]
    }

    all_data_list[[i]] <- group_data
  }

  # Combine all data using data.table's efficient rbindlist
  final_data <- rbindlist(all_data_list, use.names = TRUE, fill = TRUE)

  # Set optimal column order using data.table
  if (has_subgroups) {
    column_order <- c("simID", "group", "subgroup", "accrual.time",
                      "surv.time", "dropout.time", "tte", "total", "dropout")
  } else {
    column_order <- c("simID", "group", "accrual.time", "surv.time",
                      "dropout.time", "tte", "total", "dropout")
  }

  setcolorder(final_data, column_order)

  return(final_data)
}

# Helper function: Check if N has subgroup structure
check_subgroup_structure <- function(N) {
  if (!is.list(N)) return(FALSE)

  # Check if any element in N is a named vector (indicating subgroups)
  for (group in N) {
    if (is.vector(group) && !is.null(names(group)) && length(group) > 1 && !is.list(group)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# Helper function: Create parameter table for standard scenarios
create_parameter_table_standard <- function(N, a.time, intensity, proportion,
                                            e.time, e.hazard, d.time, d.hazard,
                                            has_subgroups, has_dropout) {

  # Extract combinations into data.table format
  combinations <- extract_combinations_standard(N, has_subgroups)

  # Create parameter data.table
  param_dt <- data.table(
    group = combinations$group,
    subgroup = combinations$subgroup,
    combination_key = combinations$combination_key
  )

  # Add parameter values using data.table operations
  param_dt[, `:=`(
    n_val = extract_parameter_values_standard(N, combination_key, has_subgroups),
    a_time_val = list(a.time),  # Common accrual time for all groups
    e_time_val = extract_parameter_values_standard(e.time, combination_key, has_subgroups),
    e_hazard_val = extract_parameter_values_standard(e.hazard, combination_key, has_subgroups)
  )]

  # Handle dropout parameters
  if (has_dropout) {
    param_dt[, `:=`(
      d_time_val = extract_parameter_values_standard(d.time, combination_key, has_subgroups),
      d_hazard_val = extract_parameter_values_standard(d.hazard, combination_key, has_subgroups)
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

# Helper function: Create parameter table for delayed subgroup scenarios
create_parameter_table_delayed <- function(N, a.time, overall.time, overall.intensity,
                                           e.time, e.hazard, d.time, d.hazard,
                                           has_subgroups, has_dropout) {

  # Extract combinations into data.table format
  combinations <- extract_combinations_standard(N, has_subgroups)

  # Create parameter data.table
  param_dt <- data.table(
    group = combinations$group,
    subgroup = combinations$subgroup,
    combination_key = combinations$combination_key
  )

  # Add parameter values using data.table operations
  param_dt[, `:=`(
    n_val = extract_parameter_values_standard(N, combination_key, has_subgroups),
    a_time_val = extract_parameter_values_delayed(a.time, combination_key, has_subgroups),
    e_time_val = extract_parameter_values_standard(e.time, combination_key, has_subgroups),
    e_hazard_val = extract_parameter_values_standard(e.hazard, combination_key, has_subgroups)
  )]

  # Handle dropout parameters
  if (has_dropout) {
    param_dt[, `:=`(
      d_time_val = extract_parameter_values_standard(d.time, combination_key, has_subgroups),
      d_hazard_val = extract_parameter_values_standard(d.hazard, combination_key, has_subgroups)
    )]
  } else {
    param_dt[, `:=`(
      d_time_val = list(NULL),
      d_hazard_val = list(NULL)
    )]
  }

  # Calculate proportions for delayed subgroup accrual
  proportions_dt <- calculate_delayed_proportions(N, a.time, overall.time, overall.intensity, has_subgroups)

  # Merge proportions back to main parameter table
  param_dt <- param_dt[proportions_dt, on = "combination_key"]

  # Set accrual parameters
  param_dt[, `:=`(
    intensity_val = list(NULL),
    proportion_val = proportion_val
  )]

  return(param_dt)
}

# Helper function: Extract group and subgroup combinations
extract_combinations_standard <- function(param_list, has_subgroups) {

  group_vals <- integer(0)
  subgroup_vals <- character(0)
  combination_keys <- character(0)

  if (has_subgroups) {
    for (g in seq_along(param_list)) {
      group_data <- param_list[[g]]
      if (is.vector(group_data) && !is.null(names(group_data)) && !is.list(group_data)) {
        # Named vector with subgroups
        subgroup_names <- names(group_data)
        for (s in seq_along(group_data)) {
          group_vals <- c(group_vals, g)
          subgroup_vals <- c(subgroup_vals, subgroup_names[s])
          combination_keys <- c(combination_keys, paste0("group", g, "_", subgroup_names[s]))
        }
      } else {
        # Single value - treat as no subgroups for this group
        group_vals <- c(group_vals, g)
        subgroup_vals <- c(subgroup_vals, NA_character_)
        combination_keys <- c(combination_keys, paste0("group", g))
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

# Helper function: Extract parameter values for specific combinations
extract_parameter_values_standard <- function(param_list, combination_keys, has_subgroups) {

  result <- vector("list", length(combination_keys))

  for (i in seq_along(combination_keys)) {
    key <- combination_keys[i]

    if (has_subgroups && grepl("_", key)) {
      # Parse group and subgroup from key
      parts <- strsplit(key, "_")[[1]]
      group_idx <- as.integer(gsub("group", "", parts[1]))
      subgroup_name <- parts[2]

      group_data <- param_list[[group_idx]]

      if (is.vector(group_data) && !is.null(names(group_data)) && !is.list(group_data)) {
        # Named vector - extract by name (for N parameter)
        if (subgroup_name %in% names(group_data)) {
          extracted_value <- group_data[subgroup_name]
          result[[i]] <- as.numeric(unname(extracted_value))
        } else {
          result[[i]] <- NA
        }
      } else if (is.list(group_data)) {
        # List structure (for e.time, e.hazard, d.time, d.hazard)
        if (subgroup_name %in% names(group_data)) {
          result[[i]] <- group_data[[subgroup_name]]
        } else {
          result[[i]] <- NA
        }
      } else {
        # Single value
        result[[i]] <- group_data
      }
    } else {
      # No subgroups or simple group key
      group_idx <- as.integer(gsub("group", "", gsub("_.*", "", key)))
      group_data <- param_list[[group_idx]]

      if (is.vector(group_data) && length(group_data) == 1 && !is.list(group_data)) {
        result[[i]] <- as.numeric(unname(group_data))
      } else if (is.list(group_data) && length(group_data) == 1) {
        result[[i]] <- group_data[[1]]
      } else {
        result[[i]] <- group_data
      }
    }
  }

  return(result)
}

# Helper function: Extract parameter values for delayed scenarios
extract_parameter_values_delayed <- function(a.time, combination_keys, has_subgroups) {

  result <- vector("list", length(combination_keys))

  for (i in seq_along(combination_keys)) {
    key <- combination_keys[i]

    if (has_subgroups && grepl("_", key)) {
      # Parse subgroup from key
      parts <- strsplit(key, "_")[[1]]
      subgroup_name <- parts[2]

      if (subgroup_name %in% names(a.time)) {
        result[[i]] <- a.time[[subgroup_name]]
      } else {
        result[[i]] <- NA
      }
    } else {
      # Handle case without subgroups - use first element
      result[[i]] <- a.time[[1]]
    }
  }

  return(result)
}

# Helper function: Calculate proportions for delayed subgroup accrual
calculate_delayed_proportions <- function(N, a.time, overall.time, overall.intensity, has_subgroups) {

  # Calculate total patients for each subgroup
  total_patients_by_subgroup <- list()
  if (has_subgroups) {
    for (group_name in names(N)) {
      if (is.vector(N[[group_name]]) && !is.null(names(N[[group_name]]))) {
        for (subgroup_name in names(N[[group_name]])) {
          if (is.null(total_patients_by_subgroup[[subgroup_name]])) {
            total_patients_by_subgroup[[subgroup_name]] <- 0
          }
          total_patients_by_subgroup[[subgroup_name]] <-
            total_patients_by_subgroup[[subgroup_name]] + N[[group_name]][subgroup_name]
        }
      }
    }
  }

  # Calculate expected recruitment based on overall.intensity and overall.time
  overall_intervals <- diff(overall.time)
  total_expected_recruitment <- sum(overall.intensity * overall_intervals)

  # Calculate proportions for each combination
  combinations <- extract_combinations_standard(N, has_subgroups)
  result_dt <- data.table(
    combination_key = combinations$combination_key,
    proportion_val = vector("list", length(combinations$combination_key))
  )

  for (i in seq_along(combinations$combination_key)) {
    key <- combinations$combination_key[i]

    if (has_subgroups && grepl("_", key)) {
      # Parse subgroup from key
      parts <- strsplit(key, "_")[[1]]
      subgroup_name <- parts[2]

      if (subgroup_name %in% names(a.time)) {
        subgroup_time <- a.time[[subgroup_name]]
        n_intervals <- length(subgroup_time) - 1

        if (n_intervals > 0) {
          # Calculate time-weighted proportions based on overlap with overall.time
          interval_proportions <- numeric(n_intervals)

          for (j in seq_len(n_intervals)) {
            interval_start <- subgroup_time[j]
            interval_end <- subgroup_time[j + 1]

            # Find overlapping periods with overall.time
            overlap_total <- 0
            for (k in seq_len(length(overall.intensity))) {
              overall_start <- overall.time[k]
              overall_end <- overall.time[k + 1]

              # Calculate overlap
              overlap_start <- max(interval_start, overall_start)
              overlap_end <- min(interval_end, overall_end)

              if (overlap_start < overlap_end) {
                overlap_duration <- overlap_end - overlap_start
                overlap_total <- overlap_total + (overall.intensity[k] * overlap_duration)
              }
            }

            interval_proportions[j] <- overlap_total
          }

          # Normalize to proportions
          if (sum(interval_proportions) > 0) {
            interval_proportions <- interval_proportions / sum(interval_proportions)
          } else {
            interval_proportions <- rep(1/n_intervals, n_intervals)
          }

          result_dt$proportion_val[[i]] <- interval_proportions
        } else {
          result_dt$proportion_val[[i]] <- NULL
        }
      } else {
        result_dt$proportion_val[[i]] <- NULL
      }
    } else {
      # Handle case without subgroups - use overall intensity pattern
      overall_intervals <- diff(overall.time)
      if (length(overall_intervals) > 0) {
        proportions <- (overall.intensity * overall_intervals) / sum(overall.intensity * overall_intervals)
        result_dt$proportion_val[[i]] <- proportions
      } else {
        result_dt$proportion_val[[i]] <- NULL
      }
    }
  }

  return(result_dt)
}
