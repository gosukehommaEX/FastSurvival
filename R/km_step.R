#' Kaplan-Meier curve drawn as steps
#'
#' Step coordinates for a Kaplan-Meier curve, including the flat segment at 1
#' and the flat tail from the last event time to the largest observed time.
#'
#' @param km The output of \code{km_one_group()} for one group.
#' @param col,lty,lwd Colour, type, and width of the line.
#' @return Called to draw; the value is not used.
#' @noRd
km_step <- function(km, col, lty, lwd) {
  xs <- step_extend(km$te, km$tmax)
  ys <- c(1, km$surv)
  if (length(xs) > length(ys)) ys <- c(ys, ys[length(ys)])
  graphics::lines(xs, ys, type = "s", col = col, lty = lty, lwd = lwd)
}
