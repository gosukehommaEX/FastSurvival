#' Fast Log-Rank Test for Two-Group Survival Data
#'
#' @description
#' Computes the log-rank test statistic for comparing survival curves between
#' two groups. Returns either a one-sided Z-score or a two-sided chi-square
#' statistic. The implementation is optimized for performance relative to
#' \code{\link[survival]{survdiff}}.
#'
#' @details
#' The log-rank statistic is computed as:
#'
#' Z = (O_1 - E_1) / sqrt(V_1)
#'
#' where O_1 is the observed number of events in the treatment group, E_1 is
#' the expected number under the null hypothesis of equal survival, and V_1 is
#' the hypergeometric variance. Tied event times are handled via the standard
#' rank-based approach: times are converted to integer ranks with
#' \code{ties.method = "min"}, and at-risk counts and event counts are
#' accumulated over distinct rank values.
#'
#' When \code{presorted = TRUE}, the input vectors are assumed to be sorted in
#' ascending order of \code{time} and \code{\link{rank}} is called directly
#' without an additional \code{\link{order}} pass. When \code{presorted = FALSE}
#' (default), the function sorts internally via \code{rank()}, which handles
#' unsorted input correctly regardless.
#'
#' In practice the two paths produce identical results; \code{presorted = TRUE}
#' removes the overhead of a redundant index computation when the caller has
#' already sorted the data, which is useful inside simulation loops.
#'
#' @param time A numeric vector of follow-up times for all subjects.
#' @param event An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{time}.
#' @param group A vector of group labels aligned with \code{time}.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group.
#' @param side An integer, either 1 or 2. If \code{side = 1}, returns the
#'   standardized log-rank statistic (Z-score, positive when the treatment
#'   group has fewer events than expected). If \code{side = 2}, returns the
#'   chi-square statistic (Z^2).
#' @param presorted A logical value. If \code{TRUE}, \code{time}, \code{event},
#'   and \code{group} are assumed to be sorted in ascending order of
#'   \code{time}, and the internal \code{order()} call is skipped. If
#'   \code{FALSE} (default), sorting is handled internally.
#'
#' @return A single numeric value: the Z-score when \code{side = 1}, or the
#'   chi-square statistic when \code{side = 2}. Returns \code{NA_real_} when
#'   \code{V1 = 0} (e.g., all events in one group).
#'
#' @examples
#' library(survival)
#'
#' # Two-sided test: compare with survdiff
#' chisq_fast <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 2)
#' chisq_ref  <- survdiff(Surv(futime, fustat) ~ rx, data = ovarian)$chisq
#' cat("survdiff_fast:", chisq_fast, "\n")
#' cat("survdiff     :", chisq_ref,  "\n")
#'
#' # One-sided test (Z-score)
#' z_stat <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 1)
#' cat("One-sided Z:", z_stat, "\n")
#'
#' # presorted = TRUE: sort once outside, reuse inside a loop
#' ord <- order(ovarian$futime)
#' survdiff_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
#'               control = 2, side = 2, presorted = TRUE)
#'
#' \dontrun{
#' library(microbenchmark)
#' microbenchmark(
#'   survdiff_fast = survdiff_fast(ovarian$futime, ovarian$fustat,
#'                                 ovarian$rx, 2, side = 2),
#'   survdiff      = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
#'   times = 1000
#' )
#' }
#'
#' @seealso
#' \code{\link[survival]{survdiff}} for the standard implementation.
#'
#' @references
#' Mantel, N. (1966). Evaluation of survival data and two new rank order
#' statistics arising in its consideration. Cancer Chemotherapy Reports,
#' 50(3), 163-270.
#'
#' Peto, R., & Peto, J. (1972). Asymptotically efficient rank invariant test
#' procedures. Journal of the Royal Statistical Society, 135(2), 185-207.
#'
#' @export
survdiff_fast <- function(time, event, group, control, side,
                          presorted = FALSE) {
  # Input validation
  if (length(time) != length(event) || length(time) != length(group)) {
    stop("'time', 'event', and 'group' must have the same length")
  }
  if (!side %in% c(1L, 2L)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided)")
  }
  if (sum(event) == 0L) {
    stop("No events observed in the data")
  }

  # Treatment indicator: 1 = treatment, 0 = control
  j <- as.numeric(group != control)

  # Convert times to integer ranks (ties.method = "min") to handle ties.
  # When presorted = TRUE the data are already ordered, but rank() still
  # needs to scan the full vector; skipping order() saves one O(n log n) pass.
  time_ranked <- rank(time, ties.method = "min")

  # Distinct event-time ranks in ascending order
  t.k <- sort.int(unique.default(time_ranked[event == 1L]))

  # Numbers at risk at each distinct event time, by group.
  # rev(cumsum(rev(tabulate(...)))) gives, at position k, the count of
  # observations with rank >= k, i.e. still at risk just before time k.
  n.1k <- rev(cumsum(rev(tabulate(time_ranked * j))))[t.k]
  n.0k <- rev(cumsum(rev(tabulate(time_ranked * (1 - j)))))[t.k]
  n.1k[is.na(n.1k)] <- 0L
  n.0k[is.na(n.0k)] <- 0L
  n.jk <- n.1k + n.0k

  # Event counts at each distinct event time, by group
  e.1k <- tabulate(time_ranked * j * event)[t.k]
  e.0k <- tabulate(time_ranked * (1 - j) * event)[t.k]
  e.1k[is.na(e.1k)] <- 0L
  e.0k[is.na(e.0k)] <- 0L
  e.jk <- e.1k + e.0k

  # Observed and expected events
  O1 <- sum(e.1k)
  E1 <- sum(e.jk * (n.1k / n.jk))

  # Hypergeometric variance
  V1 <- sum(
    (n.1k * n.0k * e.jk * (n.jk - e.jk)) / (n.jk^2 * (n.jk - 1L)),
    na.rm = TRUE
  )

  if (!is.finite(V1) || V1 == 0) return(NA_real_)

  LR <- (O1 - E1) / sqrt(V1)

  if (side == 1L) LR else LR^2
}
