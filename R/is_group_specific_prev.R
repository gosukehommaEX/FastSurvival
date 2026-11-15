#' Detect group-specific prevalence
#'
#' @param prev The \code{prevalence} argument of \code{simdata_fast()}.
#' @return \code{TRUE} when \code{prev} is a list with exactly the elements
#'   \code{control} and \code{treatment}.
#' @noRd
is_group_specific_prev <- function(prev) {
  is.list(prev) && !is.array(prev) && !is.null(names(prev)) &&
    setequal(names(prev), c("control", "treatment"))
}
