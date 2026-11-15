#' Build the one-row summary for a single scenario
#'
#' @param s One scenario of a \code{scenario_fast} object, a list with
#'   \code{args}, \code{label}, and \code{null}.
#' @param tmax End of the time grid.
#' @return A one-row data frame with the scenario, the null flag, the sample
#'   size, the two medians, the hazard ratios at the second and the last grid
#'   points, and whether the curves cross.
#' @noRd
scenario_summary_row <- function(s, tmax) {
  a     <- s$args
  times <- seq(0, tmax, length.out = 201L)
  ev    <- scenario_eval(a, times)

  n_arg <- a[["n"]]
  n_tot <- if (is.null(n_arg)) NA_integer_ else as.integer(sum(n_arg))

  med_c <- curve_median(ev$times, ev$s_c)
  med_t <- if (isTRUE(ev$two)) curve_median(ev$times, ev$s_t) else NA_real_

  if (isTRUE(ev$two)) {
    hr_start <- ev$hr[2L]
    hr_end   <- ev$hr[length(ev$hr)]
    diff_st  <- ev$s_t - ev$s_c
    crossing <- any(diff_st[-1L] > 0) && any(diff_st[-1L] < 0)
  } else {
    hr_start <- NA_real_
    hr_end   <- NA_real_
    crossing <- NA
  }

  data.frame(
    Scenario = s$label,
    Null     = s$null,
    N        = n_tot,
    Median_C = round(med_c, 2),
    Median_T = round(med_t, 2),
    HR_start = round(hr_start, 3),
    HR_end   = round(hr_end, 3),
    Crossing = crossing,
    stringsAsFactors = FALSE
  )
}
