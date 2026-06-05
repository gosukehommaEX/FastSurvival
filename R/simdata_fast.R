#' Fast Simulation of Two-Group Time-to-Event Trial Data
#'
#' @description
#' Simulates time-to-event trial data for one or two groups across many
#' simulated trials, with piecewise accrual, piecewise-exponential survival and
#' dropout, and optional subgroups defined by a prevalence specification. The
#' entire generation pipeline (accrual, survival, dropout, derived columns, and
#' two-group interleaving) runs in a single C++ kernel that materializes the
#' output data frame once, avoiding intermediate R-level vector operations and
#' copies. The random-number stream is consumed in the same order as a per-group
#' reference implementation, so results are reproducible from \code{seed}.
#'
#' @details
#' For each subject the observed time-to-event is
#' \code{tte = pmin(surv_time, dropout_time)} and \code{event} is 1 when the
#' survival time occurs first. The calendar time of the observed event is
#' \code{accrual_time + tte}.
#'
#' The total enrolled is fixed at \code{sum(n)}. With \code{a.rate} the rates are
#' absolute (subjects per unit time): when the accrual period is fully specified
#' the rates must accrue exactly \code{sum(n)}, and when one extra rate is given
#' the end of the final interval is solved so the total is met. With \code{a.prop}
#' the values are relative proportions that distribute \code{sum(n)} across the
#' fully specified intervals. Each accrual interval receives a deterministic
#' number of subjects (the rate or proportion times the group total, rounded to
#' keep the per-group total exact), placed uniformly within the interval.
#'
#' Survival and dropout are exponential when a single hazard (or median) is
#' supplied and piecewise-exponential when a vector is supplied together with
#' the corresponding \code{e.time} or \code{d.time} breakpoints, whose last
#' element must be \code{Inf}. Group-specific parameters are supplied as a
#' two-element list (control first, treatment second).
#'
#' When \code{prevalence} is supplied the trial has subgroups. A numeric vector
#' defines a single factor; a list of numeric vectors defines several
#' independent factors; a multi-dimensional array defines the joint
#' distribution of correlated factors. Per-cell hazards may be supplied as a
#' list with one element per cell. With \code{fixed.alloc = TRUE} the subgroup
#' sizes are deterministic; otherwise subgroup membership is drawn from the
#' prevalence distribution.
#'
#' @param nsim Number of simulated trials.
#' @param n Either a single total sample size (split by \code{alloc}) or a
#'   length-two vector of per-group sample sizes.
#' @param alloc A length-two allocation ratio, used when \code{n} is scalar.
#' @param a.time A numeric vector of accrual-interval breakpoints.
#' @param a.rate Absolute accrual rates (subjects per unit time), interpreted in
#'   one of two ways. With length \code{length(a.time) - 1} the accrual period is
#'   fully specified and the rates must accrue exactly \code{sum(n)} subjects (an
#'   inconsistent total is an error). With length \code{length(a.time)} the final
#'   rate applies to an open last interval whose end time is computed so the
#'   total is \code{sum(n)}. Supply exactly one of \code{a.rate} and \code{a.prop}.
#' @param a.prop Accrual proportions, one per accrual interval (length
#'   \code{length(a.time) - 1}), giving the fraction of subjects enrolled in each
#'   interval. Values are normalized to sum to one and distribute the fixed total
#'   \code{sum(n)}. Unlike \code{a.rate} this carries no rate, so the accrual
#'   period must be fully specified by \code{a.time}. Supply exactly one of
#'   \code{a.rate} and \code{a.prop}.
#' @param e.hazard Survival hazard(s). A scalar or vector for one group, or a
#'   two-element list for two groups; per-cell lists are used with subgroups.
#' @param e.median Survival median(s); an alternative to \code{e.hazard}.
#' @param e.time Survival breakpoints for piecewise hazards (last element
#'   \code{Inf}).
#' @param d.hazard Dropout hazard(s), same structure as \code{e.hazard}.
#' @param d.median Dropout median(s); an alternative to \code{d.hazard}.
#' @param d.time Dropout breakpoints for piecewise hazards.
#' @param seed Optional integer seed for the \code{dqrng} generator.
#' @param prevalence Optional subgroup prevalence specification (numeric
#'   vector, list of vectors, array, or a named \code{control}/\code{treatment}
#'   list for group-specific prevalence).
#' @param fixed.alloc Logical; when \code{TRUE} subgroup sizes are
#'   deterministic rather than drawn.
#'
#' @return A \code{data.frame} with \code{nsim * sum(n)} rows. The columns are
#'   \code{sim}, \code{group}, any subgroup columns, \code{accrual_time},
#'   \code{surv_time}, \code{dropout_time}, \code{tte}, \code{event}, and
#'   \code{calendar_time}.
#'
#' @examples
#' # One-group simulation, simple exponential, no dropout
#' df1 <- simdata_fast(
#'   nsim     = 100,
#'   n        = 50,
#'   a.time   = c(0, 12),
#'   a.rate   = 50 / 12,
#'   e.median = 18,
#'   seed     = 1
#' )
#' head(df1)
#'
#' # Accrual rate with the final interval computed from the total: 20 per unit
#' # time for the first 12 units, then 30 per unit time until 500 are enrolled
#' df1b <- simdata_fast(
#'   nsim     = 100,
#'   n        = 500,
#'   a.time   = c(0, 12),
#'   a.rate   = c(20, 30),
#'   e.median = 18,
#'   seed     = 1
#' )
#' head(df1b)
#'
#' # Accrual by proportion: 30% enrolled in [0, 6], 70% in [6, 12]
#' df1c <- simdata_fast(
#'   nsim     = 100,
#'   n        = 50,
#'   a.time   = c(0, 6, 12),
#'   a.prop   = c(0.3, 0.7),
#'   e.median = 18,
#'   seed     = 1
#' )
#' head(df1c)
#'
#' # Two-group simulation, simple exponential, with dropout
#' df2 <- simdata_fast(
#'   nsim     = 100,
#'   n        = c(60, 60),
#'   a.time   = c(0, 6, 12),
#'   a.rate   = c(8, 12),
#'   e.median = list(18, 24),
#'   d.hazard = list(0.01, 0.01),
#'   seed     = 2
#' )
#' head(df2)
#'
#' # One factor with three levels: single subgroup column
#' df3 <- simdata_fast(
#'   nsim       = 100,
#'   n          = c(150, 150),
#'   a.time     = c(0, 12),
#'   a.rate     = 300 / 12,
#'   e.hazard   = list(list(0.10, 0.08, 0.06), 0.05),
#'   prevalence = c(0.5, 0.3, 0.2),
#'   seed       = 3
#' )
#' head(df3)
#'
#' # Two independent factors (2 x 2): columns subgroup1 and subgroup2.
#' # Four cells in column-major order: (1,1), (2,1), (1,2), (2,2).
#' df4 <- simdata_fast(
#'   nsim       = 100,
#'   n          = 200,
#'   a.time     = c(0, 12),
#'   a.rate     = 200 / 12,
#'   e.hazard   = list(0.10, 0.08, 0.07, 0.05),
#'   prevalence = list(c(0.5, 0.5), c(0.6, 0.4)),
#'   seed       = 4
#' )
#' head(df4)
#'
#' # Two correlated factors via a joint-distribution array (2 x 2)
#' df5 <- simdata_fast(
#'   nsim       = 100,
#'   n          = 200,
#'   a.time     = c(0, 12),
#'   a.rate     = 200 / 12,
#'   e.hazard   = 0.08,
#'   prevalence = array(c(0.40, 0.10, 0.15, 0.35), dim = c(2, 2)),
#'   seed       = 5
#' )
#' head(df5)
#'
#' @seealso \code{\link{analysis_fast}}
#'
#' @export
simdata_fast <- function(nsim       = 1000,
                         n,
                         alloc      = c(1, 1),
                         a.time,
                         a.rate     = NULL,
                         a.prop     = NULL,
                         e.hazard   = NULL,
                         e.median   = NULL,
                         e.time     = NULL,
                         d.hazard   = NULL,
                         d.median   = NULL,
                         d.time     = NULL,
                         seed       = NULL,
                         prevalence = NULL,
                         fixed.alloc = FALSE) {

  if (!is.null(seed)) dqrng::dqset.seed(seed)

  use_subgroup        <- !is.null(prevalence)
  group_specific_prev <- use_subgroup && is_group_specific_prev(prevalence)

  spec_ctrl <- NULL
  spec_trt  <- NULL
  if (use_subgroup) {
    if (group_specific_prev) {
      spec_ctrl <- build_prevalence_spec(prevalence[["control"]])
      spec_trt  <- build_prevalence_spec(prevalence[["treatment"]])
      if (spec_ctrl$n_cell != spec_trt$n_cell ||
          !identical(dim(spec_ctrl$level_table), dim(spec_trt$level_table))) {
        stop("'control' and 'treatment' prevalence must use the same ",
             "number of factors and levels")
      }
    } else {
      spec_ctrl <- build_prevalence_spec(prevalence)
      spec_trt  <- spec_ctrl
    }
  }

  if (length(n) == 1L) {
    n_grp    <- round(n * alloc / sum(alloc))
    n_groups <- 2L
  } else if (length(n) == 2L) {
    n_grp    <- n
    n_groups <- 2L
  } else {
    stop("'n' must be a scalar (total N) or a vector of length 2 (per-group)")
  }

  if (use_subgroup) {
    is_two_group <- (length(n) == 2L) || group_specific_prev
  } else {
    is_two_group <- n_groups == 2L && (
      is.list(e.hazard) || is.list(e.median) ||
        is.list(d.hazard) || is.list(d.median)
    )
    if (!is_two_group && length(n) == 2L) {
      stop("Two-group sample sizes supplied but hazard/median parameters are not lists.\n",
           "Wrap group-specific parameters in list().")
    }
  }

  if (!is_two_group) {
    n_grp    <- if (length(n) == 1L) n else n[1L]
    n_groups <- 1L
  }

  if (use_subgroup && is_two_group) {
    chk_two <- function(x, nm) {
      if (!is.null(x) && !(is.list(x) && length(x) == 2L)) {
        stop("For a two-group simulation, '", nm,
             "' must be a list of length 2 (one element per group).")
      }
    }
    chk_two(e.hazard, "e.hazard")
    chk_two(e.median, "e.median")
    chk_two(d.hazard, "d.hazard")
    chk_two(d.median, "d.median")
  }

  if (!is.null(e.hazard) && !is.null(e.median)) {
    stop("Specify exactly one of 'e.hazard' and 'e.median', not both")
  }
  if (is.null(e.hazard) && is.null(e.median)) {
    stop("One of 'e.hazard' or 'e.median' must be supplied")
  }
  if (!is.null(e.median)) e.hazard <- convert_median_to_hazard(e.median)

  has_dropout <- !is.null(d.hazard) || !is.null(d.median)
  if (!is.null(d.hazard) && !is.null(d.median)) {
    stop("Specify exactly one of 'd.hazard' and 'd.median', not both")
  }
  if (!is.null(d.median)) d.hazard <- convert_median_to_hazard(d.median)

  use_rate <- !is.null(a.rate)
  use_prop <- !is.null(a.prop)
  if (use_rate == use_prop) {
    stop("Supply exactly one of 'a.rate' and 'a.prop'")
  }
  if (length(a.time) < 2L) stop("'a.time' must have at least two elements")
  if (any(diff(a.time) <= 0)) stop("'a.time' must be strictly increasing")

  n_total    <- sum(n_grp)
  n_int_time <- length(a.time) - 1L
  acc_tol    <- 1e-8 * max(1, n_total)

  if (use_rate) {
    if (any(a.rate <= 0)) stop("All 'a.rate' values must be positive")
    if (length(a.rate) == n_int_time) {
      implied <- sum(a.rate * diff(a.time))
      if (abs(implied - n_total) > acc_tol) {
        stop("'a.rate' implies ", round(implied, 4), " subjects over the ",
             "accrual period but sum(n) is ", n_total, ".\n",
             "Make them consistent, drop the final 'a.time' breakpoint to let ",
             "it be computed from the rate, or use 'a.prop' for relative ",
             "accrual.")
      }
      a.time_full <- a.time
    } else if (length(a.rate) == n_int_time + 1L) {
      bounded <- if (n_int_time >= 1L) {
        sum(a.rate[seq_len(n_int_time)] * diff(a.time))
      } else {
        0
      }
      remaining <- n_total - bounded
      if (remaining <= 0) {
        stop("The specified accrual intervals already accrue ",
             round(bounded, 4), " subjects, at least sum(n) = ", n_total,
             "; the final interval cannot be extended.")
      }
      final_dur   <- remaining / a.rate[length(a.rate)]
      a.time_full <- c(a.time, a.time[length(a.time)] + final_dur)
    } else {
      stop("'a.rate' must have length equal to length(a.time) - 1 (fully ",
           "specified accrual) or length(a.time) (final interval end computed ",
           "from sum(n))")
    }
    weights_a <- a.rate * diff(a.time_full)
  } else {
    if (any(a.prop <= 0)) stop("All 'a.prop' values must be positive")
    if (length(a.prop) != n_int_time) {
      stop("'a.prop' must have length equal to length(a.time) - 1")
    }
    a.time_full <- a.time
    weights_a   <- as.numeric(a.prop)
  }

  a_int_prob <- weights_a / sum(weights_a)

  n_cell_ctrl <- if (use_subgroup) spec_ctrl$n_cell else 1L
  n_cell_trt  <- if (use_subgroup) spec_trt$n_cell  else 1L
  n_cell      <- n_cell_ctrl

  # Build the per-cell hazard / piecewise pre-computation lists for one group.
  # Each returned list has one element per cell: the hazard vector, the finite
  # breakpoints, and the cumulative hazard at those breakpoints. For a single
  # hazard the breakpoint / cumulative entries are empty (unused by the kernel).
  build_exp_specs <- function(hazard_grp, time_grp, n_cell) {
    haz <- vector("list", n_cell)
    fin <- vector("list", n_cell)
    cum <- vector("list", n_cell)
    for (c in seq_len(n_cell)) {
      hz <- if (is.list(hazard_grp)) hazard_grp[[c]] else hazard_grp
      tm <- if (n_cell > 1L && is.list(time_grp)) time_grp[[c]] else time_grp
      pc <- piecewise_precompute(hz, tm)
      haz[[c]] <- pc$hazard
      fin[[c]] <- pc$fin_time
      cum[[c]] <- pc$cum_haz
    }
    list(haz = haz, fin = fin, cum = cum)
  }

  # Resolve group-level survival hazard specs for control and treatment.
  if (use_subgroup && is_two_group) {
    e_haz_c <- if (is.list(e.hazard)) e.hazard[[1L]] else e.hazard
    e_haz_t <- if (is.list(e.hazard)) e.hazard[[2L]] else e.hazard
  } else if (n_groups == 2L && !use_subgroup) {
    e_haz_c <- if (is.list(e.hazard)) e.hazard[[1L]] else e.hazard
    e_haz_t <- if (is.list(e.hazard)) e.hazard[[2L]] else e.hazard
  } else {
    e_haz_c <- e.hazard
    e_haz_t <- e.hazard
  }

  d_haz_c <- NULL; d_haz_t <- NULL
  if (has_dropout) {
    if (use_subgroup && is_two_group) {
      d_haz_c <- if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard
      d_haz_t <- if (is.list(d.hazard)) d.hazard[[2L]] else d.hazard
    } else if (n_groups == 2L && !use_subgroup) {
      d_haz_c <- if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard
      d_haz_t <- if (is.list(d.hazard)) d.hazard[[2L]] else d.hazard
    } else {
      d_haz_c <- if (use_subgroup) d.hazard else if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard
      d_haz_t <- d_haz_c
    }
  }

  e.time_c <- resolve_time_arg(e.time, 1L)
  e.time_t <- resolve_time_arg(e.time, 2L)
  d.time_c <- if (has_dropout) resolve_time_arg(d.time, 1L) else NULL
  d.time_t <- if (has_dropout) resolve_time_arg(d.time, 2L) else NULL

  e_c <- build_exp_specs(e_haz_c, e.time_c, n_cell_ctrl)
  e_t <- build_exp_specs(e_haz_t, e.time_t, n_cell_trt)
  d_c <- if (has_dropout) build_exp_specs(d_haz_c, d.time_c, n_cell_ctrl) else
    list(haz = list(numeric(0)), fin = list(numeric(0)), cum = list(numeric(0)))
  d_t <- if (has_dropout) build_exp_specs(d_haz_t, d.time_t, n_cell_trt) else
    list(haz = list(numeric(0)), fin = list(numeric(0)), cum = list(numeric(0)))

  # Subgroup descriptors.
  if (use_subgroup) {
    cum_prev_c   <- spec_ctrl$cum_prev
    cum_prev_t   <- spec_trt$cum_prev
    level_tab_c  <- matrix(as.integer(spec_ctrl$level_table),
                           nrow = nrow(spec_ctrl$level_table))
    level_tab_t  <- matrix(as.integer(spec_trt$level_table),
                           nrow = nrow(spec_trt$level_table))
    sub_names    <- spec_ctrl$sub_names
    fixed_c      <- if (fixed.alloc) fixed_cell_counts(n_grp[1L], spec_ctrl$cell_prob) else integer(n_cell_ctrl)
    fixed_t      <- if (fixed.alloc) fixed_cell_counts(if (n_groups == 2L) n_grp[2L] else n_grp[1L], spec_trt$cell_prob) else integer(n_cell_trt)
  } else {
    cum_prev_c <- numeric(0); cum_prev_t <- numeric(0)
    level_tab_c <- matrix(integer(0), nrow = 1L, ncol = 0L)
    level_tab_t <- matrix(integer(0), nrow = 1L, ncol = 0L)
    sub_names  <- character(0)
    fixed_c    <- integer(0); fixed_t <- integer(0)
  }

  n_grp_int <- if (n_groups == 1L) as.integer(n_grp[1L]) else as.integer(n_grp[1:2])

  # Deterministic per-interval accrual counts for each group. Each group's
  # counts sum to its size exactly (largest-remainder rounding), and the kernel
  # repeats them every simulation, placing subjects uniformly within intervals.
  acc_counts_c <- accrual_cell_counts(n_grp_int[1L], a_int_prob)
  acc_counts_t <- if (n_groups == 2L) {
    accrual_cell_counts(n_grp_int[2L], a_int_prob)
  } else {
    integer(0)
  }

  simdata_core_full(
    as.integer(nsim), n_grp_int,
    as.numeric(a.time_full), acc_counts_c, acc_counts_t,
    as.integer(n_cell),
    e_c$haz, e_c$fin, e_c$cum,
    e_t$haz, e_t$fin, e_t$cum,
    has_dropout,
    d_c$haz, d_c$fin, d_c$cum,
    d_t$haz, d_t$fin, d_t$cum,
    as.numeric(cum_prev_c), as.numeric(cum_prev_t),
    level_tab_c, level_tab_t,
    as.character(sub_names),
    fixed.alloc,
    as.integer(fixed_c), as.integer(fixed_t)
  )
}

# Piecewise-exponential pre-computation: validates the hazard / breakpoint
# pair and returns the hazard vector, the finite left breakpoints, and the
# cumulative hazard at those breakpoints, matching rpiece_exp_r. For a single
# hazard the breakpoint / cumulative entries are empty.
piecewise_precompute <- function(hazard, e.time) {
  if (length(hazard) == 1L) {
    return(list(hazard = as.numeric(hazard),
                fin_time = numeric(0), cum_haz = numeric(0)))
  }
  if (is.null(e.time)) {
    stop("'e.time' must be supplied when 'e.hazard' (or 'e.median') is a vector")
  }
  if (length(e.time) != length(hazard) + 1L) {
    stop("length(e.time) must equal length(e.hazard) + 1")
  }
  if (!is.infinite(e.time[length(e.time)])) {
    stop("Last element of 'e.time' must be Inf")
  }
  n_int    <- length(hazard)
  fin_time <- e.time[seq_len(n_int)]
  lengths  <- diff(e.time[seq_len(n_int)])
  cum_haz  <- c(0, cumsum(hazard[-n_int] * lengths))
  list(hazard = as.numeric(hazard),
       fin_time = as.numeric(fin_time), cum_haz = as.numeric(cum_haz))
}

# ------------------------------------------------------------------ #
#  Internal helpers: subgroup prevalence and parameter resolution
# ------------------------------------------------------------------ #

# ------------------------------------------------------------------ #
#  Internal helper: detect group-specific prevalence
# ------------------------------------------------------------------ #
is_group_specific_prev <- function(prev) {
  is.list(prev) && !is.array(prev) && !is.null(names(prev)) &&
    setequal(names(prev), c("control", "treatment"))
}

# ------------------------------------------------------------------ #
#  Internal helper: build a prevalence specification for one group
# ------------------------------------------------------------------ #
# Returns a list with the cumulative cell probabilities (cum_prev), a
# cell-by-factor level table (level_table, n_cell x n_fac), the output column
# names (sub_names), the number of cells (n_cell), and the number of factors
# (n_fac). Cells are ordered column-major (first factor varies fastest).
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

# ------------------------------------------------------------------ #
#  Internal helper: deterministic cell counts for fixed allocation
# ------------------------------------------------------------------ #
# Each cell receives floor(n * p) patients; the remaining patients (n minus
# the sum of the floors) are added one at a time starting from the first
# cell, so the last cells absorb the rounding shortfall.
fixed_cell_counts <- function(n, p) {
  base <- floor(n * p)
  rem  <- n - sum(base)
  if (rem > 0L) base[seq_len(rem)] <- base[seq_len(rem)] + 1L
  as.integer(base)
}

# ------------------------------------------------------------------ #
#  Internal helper: largest-remainder accrual counts per interval
# ------------------------------------------------------------------ #
# Distributes n subjects across accrual intervals in proportion to p, with the
# rounding remainder assigned to the intervals with the largest fractional
# parts (Hamilton's method). This keeps the per-interval counts as close as
# possible to n * p and is robust to floating-point error in p, so that e.g.
# proportions of c(0.48, 0.52) with n = 500 give exactly c(240, 260).
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

# ------------------------------------------------------------------ #
#  Internal helper: resolve a possibly group-specific time argument
# ------------------------------------------------------------------ #
resolve_time_arg <- function(time_arg, group_idx) {
  if (is.null(time_arg)) return(NULL)
  if (is.list(time_arg)) time_arg[[group_idx]] else time_arg
}

# ------------------------------------------------------------------ #
#  Internal helper: convert median survival time to hazard
# ------------------------------------------------------------------ #
convert_median_to_hazard <- function(median_arg) {
  if (is.list(median_arg)) {
    lapply(median_arg, convert_median_to_hazard)
  } else {
    log(2) / median_arg
  }
}
