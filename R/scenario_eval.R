#' Analytic survival curves and hazard ratio of a scenario
#'
#' Evaluate the analytic control and treatment survival curves and the hazard
#' ratio over a time grid for one scenario. 'two' is FALSE for a single-arm
#' spec.
#'
#' @param args The list of \code{simdata_fast()} arguments of one scenario.
#' @param times Time grid.
#' @return A list with \code{times}, the control and treatment survival
#'   \code{s_c} and \code{s_t}, the hazard ratio \code{hr} (\code{NULL} for one
#'   arm), and \code{two}.
#' @noRd
scenario_eval <- function(args, times) {
  spec_haz <- args[["e.hazard"]]
  spec_med <- args[["e.median"]]
  etime    <- args[["e.time"]]

  two <- (is.list(spec_haz) && length(spec_haz) >= 2L) ||
         (is.list(spec_med) && length(spec_med) >= 2L)

  g_c <- resolve_group_haz(spec_haz, spec_med, etime, 1L)
  s_c <- exp(-pwe_cumhaz(times, g_c$haz, g_c$etime))

  if (two) {
    g_t <- resolve_group_haz(spec_haz, spec_med, etime, 2L)
    s_t <- exp(-pwe_cumhaz(times, g_t$haz, g_t$etime))
    h_c <- pwe_hazard(times, g_c$haz, g_c$etime)
    h_t <- pwe_hazard(times, g_t$haz, g_t$etime)
    hr  <- h_t / h_c
  } else {
    s_t <- NULL
    hr  <- NULL
  }
  list(times = times, s_c = s_c, s_t = s_t, hr = hr, two = two)
}
