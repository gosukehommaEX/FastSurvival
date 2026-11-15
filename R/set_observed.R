#' Observed columns of the switched rows
#'
#' Recomputes the observed time, the event indicator, and the calendar time
#' from the latent time 't' for the rows 'rows' only, so that the observed
#' columns of the other rows are kept (for example a first endpoint censored
#' at a switch before the intermediate event by simdata_fast()). A column that
#' 'd' does not have is computed for all rows; 'cal' is NULL when no
#' calendar-time column is written. Ties count as events, and a subject with
#' neither a finite event time nor a finite dropout time is censored, as in
#' simdata_fast().
#'
#' @param d Data frame.
#' @param rows Indices of the rows to update.
#' @param t Latent event times of all rows.
#' @param drop_t Dropout times of all rows.
#' @param a Entry (accrual) times of all rows.
#' @param tte,event Names of the observed-time and event-indicator columns.
#' @param cal Name of the calendar-time column, or \code{NULL}.
#' @return \code{d} with the updated columns.
#' @noRd
set_observed <- function(d, rows, t, drop_t, a, tte, event, cal) {
  obs_t <- pmin(t, drop_t)
  obs_e <- as.integer(t < drop_t | (t == drop_t & is.finite(t)))
  if (tte %in% names(d)) {
    d[[tte]][rows] <- obs_t[rows]
  } else {
    d[[tte]] <- obs_t
  }
  if (event %in% names(d)) {
    d[[event]][rows] <- obs_e[rows]
  } else {
    d[[event]] <- obs_e
  }
  if (!is.null(cal)) {
    if (cal %in% names(d)) {
      d[[cal]][rows] <- a[rows] + obs_t[rows]
    } else {
      d[[cal]] <- a + d[[tte]]
    }
  }
  d
}
