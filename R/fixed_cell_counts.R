#' Deterministic cell counts for fixed allocation
#'
#' Each cell receives floor(n * p) patients; the remaining patients (n minus
#' the sum of the floors) are added one at a time starting from the first
#' cell, so the last cells absorb the rounding shortfall. A small tolerance in
#' the floor keeps an exact share such as 100 * 0.29 = 28.999... at 29.
#'
#' @param n Number of patients.
#' @param p Cell probabilities summing to 1.
#' @return Integer vector of the cell counts.
#' @noRd
fixed_cell_counts <- function(n, p) {
  base <- floor(n * p + 1e-8)
  rem  <- n - sum(base)
  if (rem > 0L) base[seq_len(rem)] <- base[seq_len(rem)] + 1L
  as.integer(base)
}
