#' Multi-arm (K > 2) simulation
#'
#' Generates each of the K arms with the validated single-group path of
#' simdata_fast (a scalar 'n' and a non-list hazard), sharing a common accrual
#' window supplied as fully specified proportions, and stacks the arms into one
#' data frame with a 'group' column labeled 1 to K in the order of 'n'. The
#' master seed is set by the caller and the arms consume the same dqrng stream
#' in order, so the result is reproducible from 'seed'. Subgroups and the
#' illness-death model are not supported here.
#'
#' @inheritParams simdata_fast
#' @return The data frame returned by \code{simdata_fast()} for more than two
#'   arms.
#' @noRd
simdata_fast_karm <- function(nsim, n, a.time, a.rate, a.prop,
                              e.hazard, e.median, e.time,
                              d.hazard, d.median, d.time,
                              prevalence,
                              h01.hazard, h01.median,
                              h02.hazard, h02.median) {
  K <- length(n)

  if (!is.null(prevalence)) {
    stop("Multi-arm mode (length(n) > 2) does not support 'prevalence' ",
         "subgroups.")
  }
  # The illness-death model with length(n) > 2 is rejected by simdata_fast()
  # before this function is called.
  if (!is.null(e.hazard) && !is.null(e.median)) {
    stop("Specify exactly one of 'e.hazard' and 'e.median'")
  }
  if (is.null(e.hazard) && is.null(e.median)) {
    stop("One of 'e.hazard' or 'e.median' must be supplied")
  }

  # Per-arm survival specification: a length-K list, one element per arm.
  surv_is_median <- is.null(e.hazard)
  surv_arg       <- if (surv_is_median) e.median else e.hazard
  if (!is.list(surv_arg) || length(surv_arg) != K) {
    stop("In multi-arm mode '", if (surv_is_median) "e.median" else "e.hazard",
         "' must be a list with one element per arm, so its length must equal ",
         "length(n) = ", K, ".")
  }

  # Optional per-arm dropout: NULL, a shared scalar or vector, or a length-K list.
  if (!is.null(d.hazard) && !is.null(d.median)) {
    stop("Specify at most one of 'd.hazard' and 'd.median'")
  }
  drop_is_median <- is.null(d.hazard) && !is.null(d.median)
  drop_arg       <- if (!is.null(d.hazard)) d.hazard else d.median
  has_dropout    <- !is.null(drop_arg)
  if (has_dropout && is.list(drop_arg) && length(drop_arg) != K) {
    stop("In multi-arm mode a list 'd.hazard' or 'd.median' must have one ",
         "element per arm, so its length must equal length(n) = ", K, ".")
  }

  # Resolve the whole-trial accrual once so every arm shares the same accrual
  # window and shape.
  acc         <- resolve_karm_accrual(sum(n), a.time, a.rate, a.prop)
  a.time_full <- acc$a.time_full
  a_int_prob  <- acc$a_int_prob

  # A possibly per-arm 'e.time' / 'd.time': a list is indexed per arm, anything
  # else is shared across arms.
  arm_e_time <- function(g) if (is.list(e.time)) e.time[[g]] else e.time
  arm_d_time <- function(g) if (is.list(d.time)) d.time[[g]] else d.time

  blocks <- vector("list", K)
  for (g in seq_len(K)) {
    args <- list(
      nsim   = nsim,
      n      = as.numeric(n[g]),
      a.time = a.time_full,
      a.prop = a_int_prob,
      e.time = arm_e_time(g)
    )
    if (surv_is_median) {
      args$e.median <- surv_arg[[g]]
    } else {
      args$e.hazard <- surv_arg[[g]]
    }
    if (has_dropout) {
      dv <- if (is.list(drop_arg)) drop_arg[[g]] else drop_arg
      if (drop_is_median) args$d.median <- dv else args$d.hazard <- dv
      args$d.time <- arm_d_time(g)
    }
    dfg         <- do.call(simdata_fast, args)
    dfg$group   <- g
    blocks[[g]] <- dfg
  }

  out <- do.call(rbind, blocks)
  # Group rows by simulation (then by arm) so the output matches the
  # single-endpoint contract that the analysis path relies on for its skip-sort
  # fast path.
  out <- out[order(out$sim, out$group), , drop = FALSE]
  rownames(out) <- NULL
  out
}
