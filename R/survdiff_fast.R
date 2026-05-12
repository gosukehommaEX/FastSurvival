#' Fast Log-Rank Test for Two-Group Survival Data
#'
#' @description
#' Computes the log-rank test statistic for comparing survival curves between
#' two groups. Returns either a one-sided Z-score or a two-sided chi-square
#' statistic. The C++ backend walks the pooled sorted data in a single pass,
#' maintaining per-group at-risk counters and processing tied blocks
#' atomically. The R-level group-splitting copies present in the previous
#' design are eliminated.
#'
#' @details
#' The log-rank statistic is computed as:
#'
#' Z = (O_1 - E_1) / sqrt(V_1)
#'
#' where O_1 is the observed number of events in the treatment group, E_1 is
#' the expected number under the null hypothesis of equal survival, and V_1 is
#' the hypergeometric variance. Tied event times are handled correctly: all
#' subjects sharing the same event time form a tied block, and the block is
#' processed atomically inside the C++ scan.
#'
#' When \code{presorted = TRUE}, the input vectors are assumed to be sorted in
#' ascending order of \code{time} and the internal \code{order()} call is
#' skipped. When \code{presorted = FALSE} (default), sorting is handled
#' internally. In simulation loops where the data are generated in sorted order,
#' setting \code{presorted = TRUE} avoids one O(n log n) pass.
#'
#' The C++ core (\code{logrank_core}) receives the pooled sorted data and an
#' integer group indicator directly, walks the vector once, and updates per-group
#' at-risk counts and per-tied-block statistics in a single pass. No group-split
#' vectors are allocated on the R side, eliminating four vector copies of size n
#' that the previous design required for every call.
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
#' 50(3), 163-170.
#'
#' Collett, D. (2014). Modelling Survival Data in Medical Research (3rd ed.).
#' Chapman and Hall/CRC.
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
  j <- as.integer(group != control)

  # Sort pooled data by time when not presorted
  if (!presorted) {
    ord   <- order(time)
    time  <- time[ord]
    event <- as.integer(event[ord])
    j     <- j[ord]
  } else {
    event <- as.integer(event)
  }

  # C++ core: single pass over pooled sorted data -> c(O1, E1, V1)
  ovr <- logrank_core(time, event, j)
  O1  <- ovr[1L]
  E1  <- ovr[2L]
  V1  <- ovr[3L]

  if (!is.finite(V1) || V1 == 0) return(NA_real_)

  LR <- (O1 - E1) / sqrt(V1)

  if (side == 1L) LR else LR^2
}
