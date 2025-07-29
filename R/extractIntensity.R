#' Extract accrual intensity proportions for subgroups with delayed recruitment
#'
#' This function calculates the proportion of patients recruited in each time interval
#' for each subgroup when subgroups have different recruitment periods within an
#' overall study timeline. It handles scenarios where subgroups may start recruitment
#' at different times while maintaining overall intensity targets.
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
#' The function performs the following steps:
#' \enumerate{
#'   \item Calculates total sample sizes for each subgroup across all treatment groups
#'   \item Creates unique time intervals from all accrual timelines
#'   \item Determines intensity for each interval based on overall.a.time and overall.intensity
#'   \item Identifies active subgroups for each interval
#'   \item Allocates target recruitment proportionally based on subgroup sample sizes
#'   \item Normalizes proportions so each subgroup sums to 1 across intervals
#' }
#'
#' This is particularly useful for clinical trial simulations with delayed subgroup
#' recruitment where different biomarker-defined populations may enter the study
#' at different time points while maintaining overall recruitment targets.
#'
#' @examples
#' # Example with delayed subgroup A recruitment
#' n <- list(
#'   control = c(A = 25, B = 112, C = 113),
#'   treatment = c(A = 25, B = 112, C = 113)
#' )
#' overall.a.time <- c(0, 3, 9, 12)
#' overall.intensity <- c(20, 60, 100)
#' a.time <- list(
#'   A = c(6, 9, 12),
#'   B = c(0, 3, 9, 12),
#'   C = c(0, 3, 9, 12)
#' )
#'
#' result <- extractIntensity(n, overall.a.time, overall.intensity, a.time)
#' print(result)
#'
#' # Example with three treatment groups
#' n_multi <- list(
#'   control = c(A = 20, B = 80),
#'   treatmentA = c(A = 20, B = 80),
#'   treatmentB = c(A = 20, B = 80)
#' )
#' overall.a.time_multi <- c(0, 6, 12, 18)
#' overall.intensity_multi <- c(15, 25, 20)
#' a.time_multi <- list(
#'   A = c(6, 12, 18),
#'   B = c(0, 6, 12, 18)
#' )
#'
#' result_multi <- extractIntensity(n_multi, overall.a.time_multi,
#'                                  overall.intensity_multi, a.time_multi)
#' print(result_multi)
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

  # Determine intensity for each interval based on overall.a.time
  intervals$intensity <- sapply(seq_len(nrow(intervals)), function(i) {
    interval_mid <- (intervals$start[i] + intervals$end[i]) / 2
    overall_idx <- findInterval(interval_mid, overall.a.time)
    if (overall_idx > 0 && overall_idx <= length(overall.intensity)) {
      return(overall.intensity[overall_idx])
    } else {
      return(0)
    }
  })

  # Calculate target patients per interval
  intervals$target_patients <- (intervals$end - intervals$start) * intervals$intensity

  # Function to determine active subgroups for each interval
  get_active_subgroups <- function(start_time, end_time) {
    active <- character(0)

    for (subgroup_name in names(a.time)) {
      subgroup_times <- a.time[[subgroup_name]]
      if (start_time < max(subgroup_times) && end_time > min(subgroup_times)) {
        active <- c(active, subgroup_name)
      }
    }

    return(active)
  }

  # Get subgroup names
  subgroup_names <- names(subgroup_totals)

  # Calculate proportions for each subgroup across all intervals
  # First, calculate how much each subgroup should get in each interval

  # Initialize allocation matrix: intervals x subgroups
  allocation_matrix <- matrix(0, nrow = nrow(intervals), ncol = length(subgroup_names))
  colnames(allocation_matrix) <- subgroup_names

  # Calculate total study allocation for each interval (preserving overall intensity)
  total_allocation_per_interval <- intervals$target_patients

  # For each interval, allocate to active subgroups proportionally
  for (i in seq_len(nrow(intervals))) {
    active_subgroups <- get_active_subgroups(intervals$start[i], intervals$end[i])

    if (length(active_subgroups) > 0) {
      # Get sample sizes for active subgroups
      active_n <- sapply(active_subgroups, function(sg) subgroup_totals[[sg]])
      names(active_n) <- active_subgroups

      # Allocate proportionally among active subgroups
      total_active <- sum(active_n)
      for (subgroup_name in active_subgroups) {
        allocation_matrix[i, subgroup_name] <-
          (active_n[subgroup_name] / total_active) * total_allocation_per_interval[i]
      }
    }
  }

  # Calculate total allocation for each subgroup
  subgroup_total_allocations <- colSums(allocation_matrix)

  # Convert to proportions within each subgroup (each subgroup sums to 1)
  proportion_matrix <- allocation_matrix
  for (j in seq_along(subgroup_names)) {
    if (subgroup_total_allocations[j] > 0) {
      proportion_matrix[, j] <- allocation_matrix[, j] / subgroup_total_allocations[j]
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
