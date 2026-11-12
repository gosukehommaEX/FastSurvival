#' Arguments of a Formula Call Through the Native Pipe
#'
#' With \code{d |> f(formula, ...)}, R passes the data frame \code{d} as the
#' first argument of \code{f} and the formula as the second. This internal
#' function recognizes that case, so that the estimation and testing
#' functions and \code{kmcurve_fast()} handle the call as
#' \code{f(formula, data = d, ...)}.
#'
#' @param first The first argument of the call (\code{time}, or
#'   \code{t_sorted} for \code{survfit_fast()}).
#' @param second The second argument of the call (\code{event}, or
#'   \code{e_sorted} for \code{survfit_fast()}).
#' @param data The \code{data} argument of the call.
#'
#' @return \code{NULL} when the first argument is not a data frame (or list)
#'   followed by a formula, and otherwise a list with elements \code{formula}
#'   (the second argument) and \code{data} (the first argument). Giving
#'   \code{data} as well is an error.
#'
#' @noRd
fastsurv_pipe_args <- function(first, second, data) {
  if (!is.list(first) || !inherits(second, "formula")) {
    return(NULL)
  }
  if (!is.null(data)) {
    stop("the data are given twice: through the pipe (the first argument) ",
         "and in 'data'")
  }
  list(formula = second, data = first)
}
