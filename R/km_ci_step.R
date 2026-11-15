#' Pointwise confidence limits drawn as dashed steps, of the requested type
#'
#' @param km The output of \code{km_one_group()} for one group.
#' @param conf.level Confidence level.
#' @param conf.type Type of the limits: \code{"plain"}, \code{"log"}, or
#'   \code{"log-log"}.
#' @param col Colour of the lines.
#' @return Called to draw; the value is not used.
#' @noRd
km_ci_step <- function(km, conf.level, conf.type, col) {
  if (length(km$te) == 0L) return(invisible(NULL))
  ci <- km_ci(km$surv, km$se, conf.level, conf.type)
  xs <- step_extend(km$te, km$tmax)
  lo <- c(1, ci$lower)
  up <- c(1, ci$upper)
  if (length(xs) > length(lo)) {
    lo <- c(lo, lo[length(lo)])
    up <- c(up, up[length(up)])
  }
  graphics::lines(xs, lo, type = "s", col = col, lty = 2)
  graphics::lines(xs, up, type = "s", col = col, lty = 2)
}
