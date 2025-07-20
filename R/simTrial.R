#' Simulate clinical trial data for multiple groups with optional subgroups using data.table
#'
#' This function generates simulation datasets for clinical trials with multiple treatment groups
#' and optional subgroups using data.table for maximum performance and memory efficiency.
#' It extends the functionality of simData to handle complex trial designs with time-to-event
#' endpoints, accounting for piecewise uniform distribution for patient accrual and piecewise
#' exponential distributions for patient survival time and dropout time.
#'
#' @param nsim A positive integer specifying the number of simulation iterations. Default is 1000.
#' @param N A list where each element corresponds to a group, and if subgroups exist,
#'   each element is a named vector with subgroup sample sizes. If no subgroups,
#'   a vector of sample sizes for each group.
#' @param a.time A list where each element corresponds to a group (and subgroup if applicable),
#'   containing time points defining the accrual intervals.
#' @param intensity A list where each element corresponds to a group (and subgroup if applicable),
#'   containing accrual intensities. Cannot be used with proportion. Default is NULL.
#' @param proportion A list where each element corresponds to a group (and subgroup if applicable),
#'   containing accrual proportions. Cannot be used with intensity. Default is NULL.
#' @param e.time A list where each element corresponds to a group (and subgroup if applicable),
#'   containing time points defining the survival time intervals.
#' @param e.hazard A list where each element corresponds to a group (and subgroup if applicable),
#'   containing hazard rates for survival time.
#' @param d.time A list where each element corresponds to a group (and subgroup if applicable),
#'   containing time points defining the dropout time intervals.
#' @param d.hazard A list where each element corresponds to a group (and subgroup if applicable),
#'   containing hazard rates for dropout time.
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
#' @examples
#' library(data.table)
#'
#' # Example 1: Two groups without subgroups
#' set.seed(123)
#' data1 <- simTrial(
#'   nsim = 100,
#'   N = list(group1 = 100, group2 = 100),
#'   a.time = list(
#'     group1 = c(0, 24),
#'     group2 = c(0, 24)
#'   ),
#'   intensity = list(
#'     group1 = 100/24,
#'     group2 = 100/24
#'   ),
#'   e.time = list(
#'     group1 = c(0, Inf),
#'     group2 = c(0, Inf)
#'   ),
#'   e.hazard = list(
#'     group1 = log(2) / 12,
#'     group2 = log(2) / 15
#'   ),
#'   d.time = list(
#'     group1 = c(0, Inf),
#'     group2 = c(0, Inf)
#'   ),
#'   d.hazard = list(
#'     group1 = -log(1 - 0.1) / 12,
#'     group2 = -log(1 - 0.1) / 12
#'   )
#' )
#'
#' # Example 2: Two groups with subgroups
#' data2 <- simTrial(
#'   nsim = 50,
#'   N = list(
#'     group1 = c(A = 50, B = 50),
#'     group2 = c(A = 60, B = 40)
#'   ),
#'   a.time = list(
#'     group1 = list(A = c(0, 12), B = c(0, 18)),
#'     group2 = list(A = c(0, 12), B = c(0, 18))
#'   ),
#'   intensity = list(
#'     group1 = list(A = 50/12, B = 50/18),
#'     group2 = list(A = 60/12, B = 40/18)
#'   ),
#'   e.time = list(
#'     group1 = list(A = c(0, Inf), B = c(0, Inf)),
#'     group2 = list(A = c(0, Inf), B = c(0, Inf))
#'   ),
#'   e.hazard = list(
#'     group1 = list(A = 0.05, B = 0.08),
#'     group2 = list(A = 0.03, B = 0.06)
#'   ),
#'   d.time = list(
#'     group1 = list(A = c(0, Inf), B = c(0, Inf)),
#'     group2 = list(A = c(0, Inf), B = c(0, Inf))
#'   ),
#'   d.hazard = list(
#'     group1 = list(A = 0.01, B = 0.01),
#'     group2 = list(A = 0.01, B = 0.01)
#'   )
#' )
#'
#' @import data.table
#' @export
simTrial <- function(nsim = 1e+3, N, a.time, intensity = NULL, proportion = NULL,
                     e.time, e.hazard, d.time, d.hazard, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate input parameters
  if (is.null(intensity) && is.null(proportion)) {
    stop("Either intensity or proportion must be specified")
  }

  if (!is.null(intensity) && !is.null(proportion)) {
    stop("Cannot specify both intensity and proportion arguments")
  }

  # Determine if subgroups exist
  has_subgroups <- FALSE
  if (is.list(N) && length(N) > 0) {
    first_element <- N[[1]]
    if (is.vector(first_element) && !is.null(names(first_element))) {
      has_subgroups <- TRUE
    }
  }

  # Create parameter combination table using data.table
  param_dt <- create_parameter_table(N, a.time, intensity, proportion,
                                     e.time, e.hazard, d.time, d.hazard,
                                     has_subgroups)

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
    # Ensure group_data is a data.table
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

# Internal function: Create parameter combination table using data.table
create_parameter_table <- function(N, a.time, intensity, proportion,
                                   e.time, e.hazard, d.time, d.hazard,
                                   has_subgroups) {

  # Extract combinations into data.table format
  combinations <- extract_combinations_dt(N, has_subgroups)

  # Validate parameter structure
  param_lists <- list(a.time, e.time, e.hazard, d.time, d.hazard)
  if (!is.null(intensity)) param_lists <- c(param_lists, list(intensity))
  if (!is.null(proportion)) param_lists <- c(param_lists, list(proportion))

  # Create parameter data.table
  param_dt <- data.table(
    group = combinations$group,
    subgroup = combinations$subgroup,
    combination_key = combinations$combination_key
  )

  # Add parameter values using data.table operations
  param_dt[, `:=`(
    n_val = extract_parameter_values(N, combination_key, has_subgroups),
    a_time_val = extract_parameter_values(a.time, combination_key, has_subgroups),
    e_time_val = extract_parameter_values(e.time, combination_key, has_subgroups),
    e_hazard_val = extract_parameter_values(e.hazard, combination_key, has_subgroups),
    d_time_val = extract_parameter_values(d.time, combination_key, has_subgroups),
    d_hazard_val = extract_parameter_values(d.hazard, combination_key, has_subgroups)
  )]

  if (!is.null(intensity)) {
    param_dt[, intensity_val := extract_parameter_values(intensity, combination_key, has_subgroups)]
    param_dt[, proportion_val := vector("list", nrow(param_dt))]
  } else {
    param_dt[, intensity_val := vector("list", nrow(param_dt))]
  }

  if (!is.null(proportion)) {
    param_dt[, proportion_val := extract_parameter_values(proportion, combination_key, has_subgroups)]
  } else if (!"proportion_val" %in% names(param_dt)) {
    param_dt[, proportion_val := vector("list", nrow(param_dt))]
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
