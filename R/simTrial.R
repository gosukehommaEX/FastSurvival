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
#' @param a.time A numeric vector of time points defining the accrual intervals.
#'   Time points must be in increasing order. Applied to all groups.
#' @param intensity A numeric vector of accrual intensities (rate per time unit) for each interval.
#'   Cannot be used with proportion. Applied to all groups. Default is NULL.
#' @param proportion A numeric vector of accrual proportions for each interval that sum to 1.
#'   Cannot be used with intensity. Applied to all groups. Default is NULL.
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
#'     Accrual patterns are common across all groups.
#'   }
#'   \item{Subgroup Analysis}{
#'     Within each treatment group, subgroups can be defined with distinct
#'     characteristics (e.g., different biomarker status).
#'   }
#'   \item{Common Accrual Pattern}{
#'     All groups share the same recruitment pattern as patients are randomized
#'     after enrollment.
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
#' # Example 1: Simple two-group trial without subgroups
#' set.seed(123)
#' trial1 <- simTrial(
#'   nsim = 100,
#'   N = list(group1 = 100, group2 = 100),
#'   a.time = c(0, 24),
#'   intensity = 200/24,  # Total 200 patients over 24 months
#'   e.time = list(
#'     group1 = c(0, Inf),
#'     group2 = c(0, Inf)
#'   ),
#'   e.hazard = list(
#'     group1 = log(2) / 12,    # Control: median = 12 months
#'     group2 = log(2) / 15     # Treatment: median = 15 months
#'   ),
#'   d.time = list(
#'     group1 = c(0, Inf),
#'     group2 = c(0, Inf)
#'   ),
#'   d.hazard = list(
#'     group1 = -log(1 - 0.1) / 12,  # 10% annual dropout
#'     group2 = -log(1 - 0.1) / 12   # Same dropout rate
#'   )
#' )
#'
#' # Summary of first simulation
#' print(trial1[simID == 1, .N, by = group])
#'
#' # Example 2: Trial without dropout
#' trial2 <- simTrial(
#'   nsim = 50,
#'   N = list(control = 150, treatment = 150),
#'   a.time = c(0, 6, 12, 18),
#'   intensity = c(10, 15, 20),     # Varying recruitment rates
#'   e.time = list(
#'     control = c(0, Inf),
#'     treatment = c(0, 6, Inf)     # Treatment effect after 6 months
#'   ),
#'   e.hazard = list(
#'     control = 0.1,              # Constant hazard
#'     treatment = c(0.1, 0.06)    # Reduced hazard after 6 months
#'   ),
#'   d.time = NULL,               # No dropout
#'   d.hazard = NULL,
#'   seed = 789
#' )
#'
#' # Check that dropout column is all zeros
#' print(trial2[simID == 1, .(dropout_rate = mean(dropout))])
#'
#' @seealso
#' \code{\link{simData}} for single-group simulation,
#' \code{\link{analysisData}} for creating analysis datasets,
#' \code{\link{FastLRtest}} and \code{\link{FastHRest}} for statistical analysis
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
                     e.time, e.hazard, d.time = NULL, d.hazard = NULL, seed = NULL) {

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

  # Validate dropout parameters
  if (is.null(d.time) && !is.null(d.hazard)) {
    stop("If d.time is NULL, d.hazard must also be NULL")
  }

  if (!is.null(d.time) && is.null(d.hazard)) {
    stop("If d.hazard is NULL, d.time must also be NULL")
  }

  # Check if dropout is assumed
  has_dropout <- !is.null(d.time) && !is.null(d.hazard)

  # Determine if subgroups exist
  has_subgroups <- FALSE
  if (is.list(N) && length(N) > 0) {
    first_element <- N[[1]]
    if (is.vector(first_element) && !is.null(names(first_element))) {
      has_subgroups <- TRUE
    }
  }

  # Create parameter combination table using data.table
  param_dt <- create_parameter_table_modified(N, a.time, intensity, proportion,
                                              e.time, e.hazard, d.time, d.hazard,
                                              has_subgroups, has_dropout)

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

# Internal function: Create parameter combination table using data.table (modified)
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
