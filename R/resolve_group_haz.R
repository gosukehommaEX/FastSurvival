#' Hazard and breakpoints of one group of a scenario
#'
#' Resolve one group's hazard vector and breakpoints from a scenario's args.
#' group_idx is 1 for control and 2 for treatment.
#'
#' @param spec_haz The \code{e.hazard} argument of the scenario (a vector, or a
#'   list with one element per group), or \code{NULL}.
#' @param spec_med The \code{e.median} argument of the scenario, used when
#'   \code{spec_haz} is \code{NULL}.
#' @param etime The \code{e.time} argument of the scenario (a vector, a list
#'   with one element per group, or \code{NULL}).
#' @param group_idx 1 for control and 2 for treatment.
#' @return A list with the hazard vector \code{haz} and the breakpoints
#'   \code{etime} of the group.
#' @noRd
resolve_group_haz <- function(spec_haz, spec_med, etime, group_idx) {
  if (!is.null(spec_haz)) {
    haz <- if (is.list(spec_haz)) spec_haz[[group_idx]] else spec_haz
  } else {
    med <- if (is.list(spec_med)) spec_med[[group_idx]] else spec_med
    haz <- log(2) / med
  }
  if (is.list(haz)) {
    stop("scenario plotting does not support per-subgroup survival specs")
  }
  haz <- as.numeric(haz)
  et  <- if (is.list(etime)) etime[[group_idx]] else etime
  if (length(haz) > 1L && is.null(et)) {
    stop("'e.time' is required for a piecewise (vector) hazard or median spec")
  }
  list(haz = haz, etime = et)
}
