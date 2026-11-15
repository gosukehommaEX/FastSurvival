#' Upper boundary of the area under a step curve up to tau, closed down to zero
#'
#' @param te Event times in increasing order.
#' @param surv Survival values at \code{te}.
#' @param tau Truncation time.
#' @return A list with the \code{x} and \code{y} coordinates of the polygon.
#' @noRd
step_polygon_to_tau <- function(te, surv, tau) {
  xs   <- 0
  ys   <- 1
  prev <- 1
  for (i in seq_along(te)) {
    if (te[i] > tau) break
    xs   <- c(xs, te[i], te[i])
    ys   <- c(ys, prev, surv[i])
    prev <- surv[i]
  }
  xs <- c(xs, tau)
  ys <- c(ys, prev)
  list(x = c(xs, tau, 0), y = c(ys, 0, 0))
}
