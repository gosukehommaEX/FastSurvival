#' Default time horizon of the scenario plots
#'
#' A reasonable default time horizon, derived from breakpoints, medians,
#' hazards, and accrual windows across all scenarios.
#'
#' @param scn The scenarios of a \code{scenario_fast} object, each a list with
#'   an element \code{args}.
#' @return The time horizon; 24 when no breakpoint, median, hazard, or accrual
#'   time is available.
#' @noRd
default_tmax <- function(scn) {
  cand <- numeric(0)
  for (s in scn) {
    a   <- s$args
    et  <- a[["e.time"]]
    ets <- if (is.list(et)) unlist(et) else et
    ets <- ets[is.finite(ets)]
    if (length(ets)) cand <- c(cand, max(ets) * 1.5)
    med  <- a[["e.median"]]
    meds <- if (is.list(med)) unlist(med) else med
    if (length(meds)) cand <- c(cand, max(meds) * 2)
    haz  <- a[["e.hazard"]]
    hazs <- if (is.list(haz)) unlist(haz) else haz
    hazs <- hazs[hazs > 0]
    if (length(hazs)) cand <- c(cand, (log(2) / min(hazs)) * 2)
    at  <- a[["a.time"]]
    atf <- at[is.finite(at)]
    if (length(atf)) cand <- c(cand, max(atf))
  }
  if (length(cand) == 0L) return(24)
  max(cand)
}
