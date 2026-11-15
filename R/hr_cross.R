#' First crossing of a level by a curve, by linear interpolation
#'
#' @param grid Time grid.
#' @param values Values of the curve on the grid.
#' @param level Level to be crossed.
#' @return The first crossing time, or \code{NULL} when the curve does not cross
#'   the level.
#' @noRd
hr_cross <- function(grid, values, level = 1) {
  if (length(grid) < 2L) return(NULL)
  s   <- values - level
  idx <- which(s[-1L] * s[-length(s)] < 0)
  if (length(idx) == 0L) return(NULL)
  i  <- idx[1L]
  t0 <- grid[i]
  t1 <- grid[i + 1L]
  y0 <- s[i]
  y1 <- s[i + 1L]
  t0 - y0 * (t1 - t0) / (y1 - y0)
}
