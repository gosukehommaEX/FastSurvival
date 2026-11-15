#' Resolve a possibly group-specific time argument
#'
#' @param time_arg Breakpoints: a vector, a list with one element per group, or
#'   \code{NULL}.
#' @param group_idx Index of the group.
#' @return The breakpoints of the group, or \code{NULL}.
#' @noRd
resolve_time_arg <- function(time_arg, group_idx) {
  if (is.null(time_arg)) return(NULL)
  if (is.list(time_arg)) time_arg[[group_idx]] else time_arg
}
