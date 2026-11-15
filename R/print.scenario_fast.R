#' Print Method for Scenario Sets
#'
#' @description
#' Prints a compact one-row-per-scenario summary of a \code{scenario_fast}
#' object: the scenario label, the null flag, the total sample size, the control
#' and experimental medians, the hazard ratio at the start and end of the time
#' window, and whether the two survival curves cross. The medians and hazard
#' ratios are read from the analytic survival curves, so they are defined for
#' both exponential and piecewise-exponential scenarios.
#'
#' @param x A \code{scenario_fast} object from \code{\link{gen_scenario_fast}}.
#' @param tmax Upper limit of the time window used to read medians, the end-of-
#'   window hazard ratio, and the crossing flag. Defaults to a value derived from
#'   the scenarios.
#' @param ... Further arguments, currently ignored.
#'
#' @return The input \code{x}, invisibly.
#'
#' @seealso \code{\link{gen_scenario_fast}}, \code{\link{plot.scenario_fast}}
#'
#' @export
print.scenario_fast <- function(x, tmax = NULL, ...) {
  scn   <- x$scenarios
  n_scn <- length(scn)
  if (is.null(tmax)) tmax <- default_tmax(scn)

  rows <- lapply(scn, function(s) scenario_summary_row(s, tmax))
  tab  <- do.call(rbind, rows)

  cat(sprintf("A scenario_fast object with %d scenario%s (window 0 to %s)\n\n",
              n_scn, if (n_scn == 1L) "" else "s", format(round(tmax, 1))))
  print(tab, row.names = FALSE)
  invisible(x)
}
