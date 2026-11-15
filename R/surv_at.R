#' Survival value of a step curve at arbitrary times
#'
#' @param km The output of \code{km_one_group()} for one group.
#' @param tt Times.
#' @return The survival at \code{tt}; 1 before the first event time.
#' @noRd
surv_at <- function(km, tt) {
  if (length(tt) == 0L) return(numeric(0))
  if (length(km$te) == 0L) return(rep(1, length(tt)))
  idx <- findInterval(tt, km$te)
  ifelse(idx == 0L, 1, km$surv[pmax(idx, 1L)])
}
