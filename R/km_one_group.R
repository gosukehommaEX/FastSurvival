#' Single-group Kaplan-Meier and Nelson-Aalen pieces
#'
#' Product-limit estimate with Greenwood standard error and Nelson-Aalen
#' increments for one group, computed in plain R.
#'
#' @param t_in Observed times of one group.
#' @param e_in Event indicators of one group (1 for an event, 0 for a censored
#'   time).
#' @return A list with the distinct event times \code{te}, the numbers of events
#'   \code{d} and at risk \code{Y}, the Nelson-Aalen increments \code{dLambda},
#'   the survival \code{surv} and its standard error \code{se}, the censored
#'   times \code{cens}, the sample size \code{n}, the number of events
#'   \code{n_event}, and the largest observed time \code{tmax}.
#' @noRd
km_one_group <- function(t_in, e_in) {
  n   <- length(t_in)
  ord <- order(t_in)
  tt  <- t_in[ord]
  ee  <- as.integer(e_in[ord])
  ev_t <- tt[ee == 1L]
  cens <- tt[ee == 0L]

  if (length(ev_t) == 0L) {
    return(list(te = numeric(0), d = numeric(0), Y = numeric(0),
                dLambda = numeric(0), surv = numeric(0), se = numeric(0),
                cens = cens, n = n, n_event = 0L,
                tmax = if (n > 0L) max(tt) else NA_real_))
  }

  te <- sort(unique(ev_t))
  k  <- length(te)
  d  <- as.numeric(tabulate(match(ev_t, te), nbins = k))
  Y  <- vapply(te, function(u) sum(tt >= u), numeric(1L))

  dLambda <- d / Y
  surv    <- cumprod(1 - dLambda)
  denom   <- Y * (Y - d)
  gw_term <- ifelse(denom > 0, d / denom, 0)
  se      <- surv * sqrt(cumsum(gw_term))

  list(te = te, d = d, Y = Y, dLambda = dLambda, surv = surv, se = se,
       cens = cens, n = n, n_event = length(ev_t), tmax = max(tt))
}
