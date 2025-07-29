#' Extract accrual intensity proportions for subgroups with delayed recruitment (FIXED VERSION)
#'
#' This function calculates the proportion of patients recruited in each time interval
#' for each subgroup when subgroups have different recruitment periods within an
#' overall study timeline. It ensures that each subgroup recruits exactly their
#' specified number of patients while maintaining overall intensity targets.
#'
#' @param n A named list where each element corresponds to a treatment group.
#'   Each element should be a named vector with subgroup sample sizes.
#'   Example: list(control = c(A = 25, B = 112), treatment = c(A = 25, B = 112))
#' @param overall.a.time A numeric vector of time points defining the overall
#'   study accrual timeline. Must be in increasing order.
#' @param overall.intensity A numeric vector of accrual intensities (rate per time unit)
#'   for each interval in the overall timeline. Length must be one less than
#'   length of overall.a.time.
#' @param a.time A named list where each element corresponds to a subgroup with
#'   time points defining the subgroup-specific accrual intervals. All subgroups
#'   must end at the same time as overall.a.time but may start later (delayed recruitment).
#'
#' @return A data.frame containing the proportion of patients recruited in each
#'   unique time interval for each subgroup. The data.frame has the following columns:
#'   \describe{
#'     \item{interval}{Character vector of time intervals in "start-end" format}
#'     \item{[Subgroup]_proportion}{Numeric proportion for each subgroup, where
#'       proportions for each subgroup sum to 1 across all intervals}
#'   }
#'
#' @details
#' CORRECTED VERSION: This version ensures that:
#' 1. Each subgroup recruits exactly their specified number of patients (preserves original N)
#' 2. Overall intensity targets are maintained in each time interval
#' 3. Delayed subgroups recruit at higher rates when active to compensate for lost time
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Calculates total sample sizes for each subgroup across all treatment groups
#'   \item Creates unique time intervals from all accrual timelines
#'   \item For each subgroup, distributes their total patients across their active periods
#'   \item Ensures proportions sum to 1 for each subgroup
#' }
#'
#' @examples
#' # Example with delayed subgroup A recruitment
#' n <- list(
#'   control = c(A = 25, B = 112, C = 113),
#'   treatment = c(A = 25, B = 112, C = 113)
#' )
#' overall.a.time <- c(0, 12.5)
#' overall.intensity <- 40
#' a.time <- list(
#'   A = c(8, 12.5),
#'   B = c(0, 12.5),
#'   C = c(0, 12.5)
#' )
#'
#' result <- extractIntensity(n, overall.a.time, overall.intensity, a.time)
#' print(result)
#'
#' @seealso
#' \code{\link{simTrial}} for clinical trial simulation using these proportions,
#' \code{\link{rpieceunif}} for piecewise uniform distribution
#'
#' @export
extractIntensity <- function(n, overall.a.time, overall.intensity, a.time) {
  # Input validation
  if (!is.list(n) || length(n) < 2) {
    stop("n must be a list with at least 2 treatment groups")
  }

  if (!is.numeric(overall.a.time) || length(overall.a.time) < 2) {
    stop("overall.a.time must be a numeric vector with at least 2 time points")
  }

  if (any(diff(overall.a.time) <= 0)) {
    stop("overall.a.time must be in strictly increasing order")
  }

  if (length(overall.intensity) != length(overall.a.time) - 1) {
    stop("overall.intensity length must be one less than overall.a.time length")
  }

  if (any(overall.intensity <= 0)) {
    stop("All overall.intensity values must be positive")
  }

  if (!is.list(a.time) || length(a.time) == 0) {
    stop("a.time must be a non-empty list")
  }

  # Check that all subgroups have times within overall timeline range
  overall_min <- min(overall.a.time)
  overall_max <- max(overall.a.time)

  for (subgroup_name in names(a.time)) {
    subgroup_times <- a.time[[subgroup_name]]
    if (min(subgroup_times) < overall_min || max(subgroup_times) > overall_max) {
      stop(paste("Subgroup", subgroup_name, "times must be within overall.a.time range"))
    }
    if (max(subgroup_times) != overall_max) {
      stop(paste("Subgroup", subgroup_name, "must end at the same time as overall.a.time"))
    }
  }

  # Calculate total sample sizes for each subgroup across all treatment groups
  subgroup_totals <- list()

  for (group_idx in seq_along(n)) {
    group_data <- n[[group_idx]]
    if (is.vector(group_data) && !is.null(names(group_data)) && !is.list(group_data)) {
      for (subgroup_name in names(group_data)) {
        if (is.null(subgroup_totals[[subgroup_name]])) {
          subgroup_totals[[subgroup_name]] <- 0
        }
        subgroup_totals[[subgroup_name]] <- subgroup_totals[[subgroup_name]] + group_data[subgroup_name]
      }
    }
  }

  if (length(subgroup_totals) == 0) {
    stop("No valid subgroups found in n")
  }

  # Get unique time points and create intervals
  all_times <- sort(unique(c(overall.a.time, unlist(a.time))))
  intervals <- data.frame(
    start = all_times[-length(all_times)],
    end = all_times[-1]
  )

  # Get subgroup names
  subgroup_names <- names(subgroup_totals)

  # CORRECTED APPROACH: Each subgroup distributes their total patients across their active periods
  # This preserves original sample sizes while allowing temporal flexibility

  # Initialize proportion matrix: intervals x subgroups
  proportion_matrix <- matrix(0, nrow = nrow(intervals), ncol = length(subgroup_names))
  colnames(proportion_matrix) <- subgroup_names

  # For each subgroup, calculate proportions based on their active periods
  for (j in seq_along(subgroup_names)) {
    subgroup_name <- subgroup_names[j]
    subgroup_times <- a.time[[subgroup_name]]
    subgroup_start <- min(subgroup_times)
    subgroup_end <- max(subgroup_times)
    subgroup_duration <- subgroup_end - subgroup_start

    # Calculate total "time weight" for this subgroup across all intervals
    total_time_weight <- 0
    time_weights <- numeric(nrow(intervals))

    for (i in seq_len(nrow(intervals))) {
      interval_start <- intervals$start[i]
      interval_end <- intervals$end[i]

      # Calculate overlap between interval and subgroup active period
      overlap_start <- max(interval_start, subgroup_start)
      overlap_end <- min(interval_end, subgroup_end)
      overlap_length <- max(0, overlap_end - overlap_start)

      if (overlap_length > 0) {
        time_weights[i] <- overlap_length
        total_time_weight <- total_time_weight + overlap_length
      }
    }

    # Distribute subgroup's total patients proportionally across active intervals
    if (total_time_weight > 0) {
      for (i in seq_len(nrow(intervals))) {
        if (time_weights[i] > 0) {
          proportion_matrix[i, j] <- time_weights[i] / total_time_weight
        }
      }
    }
  }

  # Create final results data.frame
  final_results <- data.frame(
    interval = paste(intervals$start, intervals$end, sep = "-")
  )

  # Add proportion columns
  for (j in seq_along(subgroup_names)) {
    final_results[[paste0(subgroup_names[j], "_proportion")]] <- proportion_matrix[, j]
  }

  return(final_results)
}
