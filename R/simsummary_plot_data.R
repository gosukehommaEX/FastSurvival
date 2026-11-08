#' Values drawn by the plot method of simsummary_fast results
#'
#' Internal helper of \code{plot.simsummary_fast()}. It checks that a result
#' of \code{simsummary_fast()} can be plotted and returns, for every block
#' (arm and population), the abscissae of the looks, the cumulative
#' probabilities of stopping for efficacy and for futility, and the
#' probabilities of stopping at each look, together with the default axis
#' limits. The plot method only draws these values, so that the values can be
#' tested.
#'
#' @param x An object of class \code{"simsummary_fast"}.
#' @param xaxis The horizontal axis: \code{"look"}, \code{"events"}, or
#'   \code{"time"} (see \code{plot.simsummary_fast()}).
#'
#' @return A list with the elements \code{blocks}, \code{looks},
#'   \code{has.futility}, \code{xcol}, \code{xlim}, and \code{ylim}.
#'   \code{blocks} has one list per block with the elements \code{label},
#'   \code{x} (the abscissae of the looks), \code{cum.efficacy} (the
#'   cumulative probability of stopping for efficacy, \code{cum.reject}),
#'   \code{cum.futility} (the cumulative sum of \code{prob.stop.futility},
#'   or \code{NULL} without futility boundaries), \code{stage} (a matrix with
#'   one column per look and the row \code{efficacy}, followed, with futility
#'   boundaries, by the row \code{futility}), and \code{names} (the labels of
#'   the bars). \code{looks} are the labels of the looks, \code{xcol} the
#'   column that gives the abscissae (\code{NULL} for \code{xaxis = "look"}),
#'   \code{xlim} the horizontal limits of the cumulative plot, and
#'   \code{ylim} a list with the default vertical limits \code{cumulative}
#'   and \code{stage}.
#'
#' @keywords internal
#' @noRd
simsummary_plot_data <- function(x, xaxis = c("look", "events", "time")) {
  xaxis <- match.arg(xaxis)

  # A result of FastSurvival 1.1.0 stores the boundary settings without the
  # looks and rows that the blocks are checked against
  bd <- attr(x, "boundary")
  if (!is.null(bd) && is.null(bd$looks)) {
    stop("this result was created by an earlier version of FastSurvival; ",
         "call simsummary_fast() again to plot it")
  }
  blocks <- simsummary_blocks(x)
  if (is.null(blocks)) {
    stop("plot() needs a simsummary_fast() result whose blocks keep all ",
         "their looks in the original order")
  }
  has_fut <- identical(bd$mode, "z") && !is.null(bd$futility)

  xcol <- switch(xaxis, look = NULL, events = "n.event.mean",
                 time = "cutoff.mean")
  if (!is.null(xcol) && !(xcol %in% names(blocks[[1L]]))) {
    stop("xaxis = \"", xaxis, "\" needs the column '", xcol, "', which the ",
         "result does not have")
  }
  if (!is.null(xcol) &&
      !any(is.finite(unlist(lapply(blocks, function(b) b[[xcol]]))))) {
    stop("xaxis = \"", xaxis, "\" needs finite values in the column '", xcol,
         "'")
  }

  # Abscissae of the looks of a block
  x_of <- function(b) {
    if (is.null(xcol)) seq_len(nrow(b)) else b[[xcol]]
  }
  # Default upper limit of the vertical axis: slightly above the largest
  # value drawn, and at most 1
  top_of <- function(v) {
    m <- suppressWarnings(max(unlist(v), na.rm = TRUE))
    if (!is.finite(m) || m <= 0) 1 else min(1, 1.08 * m)
  }

  out <- lapply(blocks, function(b) {
    stage <- if (has_fut) {
      rbind(efficacy = b$prob.stop.efficacy, futility = b$prob.stop.futility)
    } else {
      rbind(efficacy = b$prob.stop.efficacy)
    }
    list(label        = attr(b, "label"),
         x            = x_of(b),
         cum.efficacy = b$cum.reject,
         cum.futility = if (has_fut) cumsum(b$prob.stop.futility) else NULL,
         stage        = stage,
         names        = if (is.null(xcol)) {
           bd$looks
         } else {
           format(round(b[[xcol]], 1), trim = TRUE)
         })
  })

  xs   <- unlist(lapply(out, function(b) b$x))
  xlim <- range(xs[is.finite(xs)])
  if (!(diff(xlim) > 0)) xlim <- xlim + c(-0.5, 0.5)

  cum_drawn <- lapply(out, function(b) c(b$cum.efficacy, b$cum.futility))
  stage_sum <- lapply(blocks, function(b) {
    b$prob.stop.efficacy + b$prob.stop.futility
  })

  list(blocks       = out,
       looks        = bd$looks,
       has.futility = has_fut,
       xcol         = xcol,
       xlim         = xlim,
       ylim         = list(cumulative = c(0, top_of(cum_drawn)),
                           stage      = c(0, top_of(stage_sum))))
}
