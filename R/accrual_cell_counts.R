#' Largest-remainder accrual counts per interval
#'
#' Distributes n subjects across accrual intervals in proportion to p, with the
#' rounding remainder assigned to the intervals with the largest fractional
#' parts (Hamilton's method). This keeps the per-interval counts as close as
#' possible to n * p and is robust to floating-point error in p, so that e.g.
#' proportions of c(0.48, 0.52) with n = 500 give exactly c(240, 260).
#'
#' @param n Number of subjects.
#' @param p Proportions of the intervals summing to 1.
#' @return Integer vector of the counts per interval.
#' @noRd
accrual_cell_counts <- function(n, p) {
  target <- n * p
  base   <- floor(target)
  rem    <- n - sum(base)
  if (rem > 0L) {
    take <- order(target - base, decreasing = TRUE)[seq_len(rem)]
    base[take] <- base[take] + 1L
  }
  as.integer(base)
}
