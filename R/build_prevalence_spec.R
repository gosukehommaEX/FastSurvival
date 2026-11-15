#' Build a prevalence specification for one group
#'
#' Returns a list with the cumulative cell probabilities (cum_prev), a
#' cell-by-factor level table (level_table, n_cell x n_fac), the output column
#' names (sub_names), the number of cells (n_cell), and the number of factors
#' (n_fac). Cells are ordered column-major (first factor varies fastest).
#'
#' @param prev Prevalence of one group: a numeric vector (one factor), a list of
#'   numeric vectors (independent factors), or an array (joint distribution of
#'   the factors).
#' @return The list described above; it also holds the cell probabilities
#'   \code{cell_prob}.
#' @noRd
build_prevalence_spec <- function(prev) {
  if (is.array(prev) && length(dim(prev)) >= 2L) {
    # Joint distribution of correlated factors
    dims <- dim(prev)
    cell_prob <- as.numeric(prev)
    if (any(cell_prob <= 0)) stop("All 'prevalence' values must be positive")
    cell_prob   <- cell_prob / sum(cell_prob)
    level_table <- as.matrix(expand.grid(lapply(dims, seq_len)))
    n_fac       <- length(dims)
  } else if (is.list(prev)) {
    # Several independent factors (marginal prevalence per factor)
    bad <- vapply(prev, function(x) !is.numeric(x) || is.list(x), logical(1L))
    if (any(bad)) {
      stop("Each factor in a 'prevalence' list must be a numeric vector")
    }
    if (any(unlist(prev) <= 0)) stop("All 'prevalence' values must be positive")
    margins     <- lapply(prev, function(x) x / sum(x))
    levels_f    <- vapply(margins, length, integer(1L))
    level_table <- as.matrix(expand.grid(lapply(levels_f, seq_len)))
    cell_prob   <- apply(level_table, 1L, function(idx) {
      prod(mapply(function(m, i) m[i], margins, idx))
    })
    n_fac <- length(margins)
  } else {
    # Single factor
    if (!is.numeric(prev)) stop("'prevalence' must be numeric")
    if (any(prev <= 0)) stop("All 'prevalence' values must be positive")
    cell_prob   <- prev / sum(prev)
    level_table <- matrix(seq_along(prev), ncol = 1L)
    n_fac       <- 1L
  }

  n_cell    <- length(cell_prob)
  sub_names <- if (n_fac == 1L) "subgroup" else paste0("subgroup", seq_len(n_fac))
  dimnames(level_table) <- NULL

  list(cum_prev = cumsum(cell_prob), cell_prob = cell_prob,
       level_table = level_table, sub_names = sub_names,
       n_cell = n_cell, n_fac = n_fac)
}
