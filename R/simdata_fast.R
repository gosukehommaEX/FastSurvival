#' Fast Survival Data Simulation for Clinical Trials
#'
#' @description
#' Simulates individual patient data for one- or two-group time-to-event
#' trials, with optional patient subgroups defined by one or more factors.
#' Patient accrual times are drawn from a piecewise uniform distribution.
#' Survival and dropout times are drawn from either a simple exponential or a
#' piecewise exponential distribution, depending on whether a scalar or vector
#' hazard is supplied. When subgroups are requested, each patient is assigned
#' to a cell of the cross-classification of the factors by a categorical draw,
#' and survival and dropout hazards may differ across cells. Random number
#' generation uses \code{dqrng} for speed. C++ backends handle piecewise
#' sampling, cell assignment, and two-group interleaving to minimize R-level
#' overhead.
#'
#' @details
#' For each patient, three times are generated independently:
#'
#' \describe{
#'   \item{Accrual time}{Drawn from a piecewise uniform distribution defined
#'     by \code{a.time} and \code{a.rate}.}
#'   \item{Survival time}{Drawn from a simple exponential distribution when
#'     the hazard spec is a scalar, or from a piecewise exponential
#'     distribution when it is a vector (with matching \code{e.time}).}
#'   \item{Dropout time}{Drawn from the same family, governed by
#'     \code{d.hazard} / \code{d.median} and \code{d.time}, or set to
#'     \code{Inf} when dropout is not modeled.}
#' }
#'
#' The observed time-to-event is \code{tte = pmin(surv_time, dropout_time)},
#' the event indicator is \code{1} when \code{surv_time <= dropout_time}, and
#' the calendar time is \code{accrual_time + tte}.
#'
#' Exactly one of \code{e.hazard} and \code{e.median} must be supplied.
#' Exactly one of \code{d.hazard} and \code{d.median} must be supplied when
#' dropout is modeled; both must be \code{NULL} to suppress dropout.
#'
#' For a two-group trial, group-specific parameters are passed as a list of
#' length 2, where element 1 is the control group and element 2 the treatment
#' group. For a one-group trial, plain vectors or scalars are passed directly.
#'
#' Subgroups are activated by \code{prevalence}, which accepts four forms:
#'
#' \describe{
#'   \item{A numeric vector}{One factor with as many levels as the length of
#'     the vector. The output gains a single integer column \code{subgroup}.}
#'   \item{A list of numeric vectors}{Several independent factors, one vector
#'     of marginal prevalence per factor. The factors are assigned
#'     independently (their joint distribution is the product of the
#'     marginals). The output gains one integer column per factor,
#'     \code{subgroup1}, \code{subgroup2}, and so on.}
#'   \item{An array}{Several factors with a joint distribution given directly
#'     by the array, allowing the factors to be correlated. The dimensions of
#'     the array are the numbers of factor levels. The output gains one
#'     integer column per factor, \code{subgroup1}, \code{subgroup2}, and so
#'     on.}
#'   \item{A named list \code{list(control = ..., treatment = ...)}}{Two
#'     group-specific specifications, each itself one of the three forms
#'     above, so that the subgroup composition may differ between groups
#'     (as in stratified randomization). Both groups must use the same number
#'     of factors and the same number of levels per factor.}
#' }
#'
#' Cells of the factor cross-classification are ordered column-major: the
#' first factor varies fastest, matching the ordering of \code{array} and
#' \code{expand.grid} in R. Cell-specific event hazards are supplied in this
#' order. For a one-group trial, \code{e.hazard} is either a single spec
#' (broadcast to all cells) or a list with one spec per cell. For a two-group
#' trial, \code{e.hazard} is a list of length 2, and each group element is
#' itself either a single spec (broadcast to that group's cells) or a list
#' with one spec per cell. The same broadcasting rule applies to
#' \code{e.median}, \code{d.hazard}, \code{d.median}, \code{e.time}, and
#' \code{d.time}. In subgroup mode a two-group trial must be declared through
#' a length-2 \code{n} or a group-specific \code{prevalence}, so a top-level
#' hazard list is never ambiguous between groups and cells.
#'
#' @param nsim A positive integer, the number of simulation iterations.
#'   Default is 1000.
#' @param n A positive integer (one group) or a numeric vector of length 2
#'   (two groups). Alternatively a scalar total sample size together with
#'   \code{alloc}.
#' @param alloc A numeric vector of length 2 giving the allocation ratio
#'   \code{c(control, treatment)}. Used only when \code{n} is a scalar.
#' @param a.time A numeric vector of accrual interval boundaries
#'   \code{c(0, t1, ..., T)}, strictly increasing.
#' @param a.rate A numeric vector of accrual rates, of length
#'   \code{length(a.time) - 1}, all positive.
#' @param e.hazard Event hazard spec. A scalar or vector for one group, a list
#'   of two for two groups; with subgroups, a group element may be a list with
#'   one spec per cell. Cannot be used with \code{e.median}.
#' @param e.median Event median spec, same structure as \code{e.hazard},
#'   converted via \code{log(2) / median}. Cannot be used with \code{e.hazard}.
#' @param e.time Piecewise interval boundaries for the event time,
#'   \code{c(0, t1, ..., Inf)}. Same structure rules as \code{e.hazard}.
#' @param d.hazard Dropout hazard spec, same structure as \code{e.hazard}. Set
#'   \code{d.hazard} and \code{d.median} both to \code{NULL} to suppress
#'   dropout. Cannot be used with \code{d.median}.
#' @param d.median Dropout median spec. Cannot be used with \code{d.hazard}.
#' @param d.time Piecewise interval boundaries for the dropout time.
#' @param seed A single integer for reproducibility, passed to
#'   \code{dqrng::dqset.seed}. Default \code{NULL}.
#' @param prevalence Subgroup prevalence specification. \code{NULL} (default)
#'   means no subgroups. Otherwise a numeric vector (one factor), a list of
#'   numeric vectors (independent factors), an array (joint distribution of
#'   correlated factors), or a named list
#'   \code{list(control = ..., treatment = ...)} for group-specific
#'   composition. Values must be positive and are normalized internally. See
#'   Details.
#' @param fixed.alloc A logical value controlling how patients are assigned to
#'   subgroup cells. When \code{FALSE} (default), each patient is assigned at
#'   random by a categorical draw from \code{prevalence}, so cell sizes are
#'   random with the given expected proportions. When \code{TRUE}, cell sizes
#'   are fixed deterministically: each cell receives \code{floor(n * p)}
#'   patients and any remainder is distributed one per cell from the first
#'   cell onward, so the last cells absorb the rounding shortfall. For
#'   example, \code{n = 99} with \code{prevalence = c(0.5, 0.5)} gives 50 and
#'   49, and \code{n = 100} with \code{prevalence = c(1/3, 1/3, 1/3)} gives
#'   34, 33, and 33. Fixed allocation consumes no random numbers for the
#'   assignment. Ignored when \code{prevalence} is \code{NULL}.
#'
#' @return A \code{data.frame} with \code{nsim * sum(n)} rows. Without
#'   subgroups the columns are \code{sim}, \code{group}, \code{accrual_time},
#'   \code{surv_time}, \code{dropout_time}, \code{tte}, \code{event}, and
#'   \code{calendar_time}. With a single-factor \code{prevalence}, an integer
#'   column \code{subgroup} is inserted after \code{group}. With several
#'   factors, the columns \code{subgroup1}, \code{subgroup2}, and so on are
#'   inserted instead, one per factor.
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
#' # Two-group simulation, simple exponential, with dropout
#' df2 <- simdata_fast(
#'   nsim     = 100,
#'   n        = c(100, 100),
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
#' # Group-specific composition (stratified randomization)
#' df6 <- simdata_fast(
#'   nsim       = 100,
#'   n          = c(150, 150),
#'   a.time     = c(0, 12),
#'   a.rate     = 300 / 12,
#'   e.hazard   = list(0.10, 0.05),
#'   prevalence = list(control = c(0.7, 0.3), treatment = c(0.5, 0.5)),
#'   seed       = 6
#' )
#' head(df6)
#'
#' @importFrom dqrng dqset.seed dqrexp
#' @importFrom stats setNames
#' @export
simdata_fast <- function(nsim       = 1000,
                         n,
                         alloc      = c(1, 1),
                         a.time,
                         a.rate,
                         e.hazard   = NULL,
                         e.median   = NULL,
                         e.time     = NULL,
                         d.hazard   = NULL,
                         d.median   = NULL,
                         d.time     = NULL,
                         seed       = NULL,
                         prevalence = NULL,
                         fixed.alloc = FALSE) {

  # ------------------------------------------------------------------ #
  #  Seed
  # ------------------------------------------------------------------ #
  if (!is.null(seed)) dqrng::dqset.seed(seed)

  # ------------------------------------------------------------------ #
  #  Subgroup setup
  # ------------------------------------------------------------------ #
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

  # ------------------------------------------------------------------ #
  #  Resolve group sizes
  # ------------------------------------------------------------------ #
  if (length(n) == 1L) {
    n_grp    <- round(n * alloc / sum(alloc))
    n_groups <- 2L
  } else if (length(n) == 2L) {
    n_grp    <- n
    n_groups <- 2L
  } else {
    stop("'n' must be a scalar (total N) or a vector of length 2 (per-group)")
  }

  # ------------------------------------------------------------------ #
  #  Resolve group count
  # ------------------------------------------------------------------ #
  if (use_subgroup) {
    # In subgroup mode, two groups must be declared via a length-2 'n' or a
    # group-specific (named) prevalence. A top-level hazard list is otherwise
    # reserved for per-cell specs in a one-group trial.
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

  # In subgroup two-group mode require length-2 lists for group-specific
  # parameters, matching the convention used without subgroups.
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

  # ------------------------------------------------------------------ #
  #  Validate and resolve e.hazard / e.median
  # ------------------------------------------------------------------ #
  if (!is.null(e.hazard) && !is.null(e.median)) {
    stop("Specify exactly one of 'e.hazard' and 'e.median', not both")
  }
  if (is.null(e.hazard) && is.null(e.median)) {
    stop("One of 'e.hazard' or 'e.median' must be supplied")
  }
  if (!is.null(e.median)) e.hazard <- convert_median_to_hazard(e.median)

  # ------------------------------------------------------------------ #
  #  Validate and resolve d.hazard / d.median
  # ------------------------------------------------------------------ #
  has_dropout <- !is.null(d.hazard) || !is.null(d.median)
  if (!is.null(d.hazard) && !is.null(d.median)) {
    stop("Specify exactly one of 'd.hazard' and 'd.median', not both")
  }
  if (!is.null(d.median)) d.hazard <- convert_median_to_hazard(d.median)

  # ------------------------------------------------------------------ #
  #  Validate accrual parameters
  # ------------------------------------------------------------------ #
  if (length(a.rate) != length(a.time) - 1L) {
    stop("'a.rate' must have length equal to length(a.time) - 1")
  }
  if (any(a.rate <= 0)) stop("All 'a.rate' values must be positive")
  if (any(diff(a.time) <= 0)) stop("'a.time' must be strictly increasing")

  # Pre-compute accrual cumulative probability table (shared across groups)
  lengths_a <- diff(a.time)
  weights_a <- a.rate * lengths_a
  cum_p_a   <- cumsum(weights_a / sum(weights_a))

  # Subgroup cell counts (1 when no subgroups)
  n_cell_ctrl <- if (use_subgroup) spec_ctrl$n_cell else 1L
  n_cell_trt  <- if (use_subgroup) spec_trt$n_cell  else 1L

  # ------------------------------------------------------------------ #
  #  Simulate each group and combine
  # ------------------------------------------------------------------ #
  if (n_groups == 1L) {
    e.hazard_grp <- if (use_subgroup) e.hazard else if (is.list(e.hazard)) e.hazard[[1L]] else e.hazard
    d.hazard_grp <- if (!has_dropout) NULL else if (use_subgroup) d.hazard else if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard

    df <- simulate_group(
      nsim        = nsim,
      n           = n_grp,
      group_id    = 1L,
      a.time      = a.time,
      cum_p_a     = cum_p_a,
      e.hazard    = e.hazard_grp,
      e.time      = e.time,
      d.hazard    = d.hazard_grp,
      d.time      = d.time,
      n_sub       = n_cell_ctrl,
      cum_prev    = if (use_subgroup) spec_ctrl$cum_prev else NULL,
      level_table = if (use_subgroup) spec_ctrl$level_table else NULL,
      sub_names   = if (use_subgroup) spec_ctrl$sub_names else NULL,
      cell_prob   = if (use_subgroup) spec_ctrl$cell_prob else NULL,
      fixed.alloc = fixed.alloc
    )
  } else {
    e.time_ctrl <- resolve_time_arg(e.time, 1L)
    e.time_trt  <- resolve_time_arg(e.time, 2L)
    d.time_ctrl <- if (has_dropout) resolve_time_arg(d.time, 1L) else NULL
    d.time_trt  <- if (has_dropout) resolve_time_arg(d.time, 2L) else NULL

    df_ctrl <- simulate_group(
      nsim        = nsim,
      n           = n_grp[1L],
      group_id    = 1L,
      a.time      = a.time,
      cum_p_a     = cum_p_a,
      e.hazard    = if (is.list(e.hazard)) e.hazard[[1L]] else e.hazard,
      e.time      = e.time_ctrl,
      d.hazard    = if (has_dropout) (if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard) else NULL,
      d.time      = d.time_ctrl,
      n_sub       = n_cell_ctrl,
      cum_prev    = if (use_subgroup) spec_ctrl$cum_prev else NULL,
      level_table = if (use_subgroup) spec_ctrl$level_table else NULL,
      sub_names   = if (use_subgroup) spec_ctrl$sub_names else NULL,
      cell_prob   = if (use_subgroup) spec_ctrl$cell_prob else NULL,
      fixed.alloc = fixed.alloc
    )

    df_trt <- simulate_group(
      nsim        = nsim,
      n           = n_grp[2L],
      group_id    = 2L,
      a.time      = a.time,
      cum_p_a     = cum_p_a,
      e.hazard    = if (is.list(e.hazard)) e.hazard[[2L]] else e.hazard,
      e.time      = e.time_trt,
      d.hazard    = if (has_dropout) (if (is.list(d.hazard)) d.hazard[[2L]] else d.hazard) else NULL,
      d.time      = d.time_trt,
      n_sub       = n_cell_trt,
      cum_prev    = if (use_subgroup) spec_trt$cum_prev else NULL,
      level_table = if (use_subgroup) spec_trt$level_table else NULL,
      sub_names   = if (use_subgroup) spec_trt$sub_names else NULL,
      cell_prob   = if (use_subgroup) spec_trt$cell_prob else NULL,
      fixed.alloc = fixed.alloc
    )

    if (use_subgroup && n_cell_ctrl > 1L) {
      # Subgroup output carries extra columns, so interleave in R.
      df <- rbind(df_ctrl, df_trt)
      df <- df[order(df$sim, df$group), , drop = FALSE]
      rownames(df) <- NULL
    } else {
      # C++ interleave: replaces rbind() + order()
      df <- interleave_groups(df_ctrl, df_trt, nsim,
                              as.integer(n_grp[1L]), as.integer(n_grp[2L]))
    }
  }

  df
}

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
#  Internal helper: simulate one group
# ------------------------------------------------------------------ #
simulate_group <- function(nsim, n, group_id, a.time, cum_p_a,
                           e.hazard, e.time, d.hazard, d.time,
                           n_sub = 1L, cum_prev = NULL,
                           level_table = NULL, sub_names = NULL,
                           cell_prob = NULL, fixed.alloc = FALSE) {
  total_n <- nsim * n

  # Accrual times via C++ piecewise uniform sampler
  accrual_time <- rpiece_unif_cpp(total_n, a.time, cum_p_a)

  if (n_sub == 1L) {
    # No subgroups: original code path (unchanged output and RNG order)
    surv_time <- rpiece_exp_r(total_n, e.hazard, e.time)

    dropout_time <- if (!is.null(d.hazard)) {
      rpiece_exp_r(total_n, d.hazard, d.time)
    } else {
      rep(Inf, total_n)
    }

    tte           <- pmin(surv_time, dropout_time)
    event         <- as.integer(surv_time <= dropout_time)
    calendar_time <- accrual_time + tte

    return(data.frame(
      sim           = rep(seq_len(nsim), each = n),
      group         = group_id,
      accrual_time  = accrual_time,
      surv_time     = surv_time,
      dropout_time  = dropout_time,
      tte           = tte,
      event         = event,
      calendar_time = calendar_time
    ))
  }

  # Validate cell-specific spec lengths
  if (is.list(e.hazard) && length(e.hazard) != n_sub) {
    stop("Cell-specific 'e.hazard' list must have one element per cell (",
         n_sub, ")")
  }
  if (is.list(d.hazard) && length(d.hazard) != n_sub) {
    stop("Cell-specific 'd.hazard' list must have one element per cell (",
         n_sub, ")")
  }

  # Subgroup cell path: assign cells, then draw survival/dropout per cell.
  # Random assignment uses a categorical draw; fixed assignment uses
  # deterministic per-cell counts (floor with a front-loaded remainder) and
  # consumes no random numbers.
  if (fixed.alloc) {
    counts   <- fixed_cell_counts(n, cell_prob)
    cell_one <- rep.int(seq_len(n_sub), counts)
    cell     <- rep(cell_one, times = nsim)
  } else {
    cell <- rcat_cpp(total_n, cum_prev)
  }
  surv_time <- numeric(total_n)
  dropout_time <- rep(Inf, total_n)

  get_sub_spec <- function(spec, s) if (is.list(spec)) spec[[s]] else spec

  for (s in seq_len(n_sub)) {
    idx <- which(cell == s)
    if (length(idx) == 0L) next
    surv_time[idx] <- rpiece_exp_r(length(idx),
                                   get_sub_spec(e.hazard, s),
                                   get_sub_spec(e.time, s))
    if (!is.null(d.hazard)) {
      dropout_time[idx] <- rpiece_exp_r(length(idx),
                                        get_sub_spec(d.hazard, s),
                                        get_sub_spec(d.time, s))
    }
  }

  tte           <- pmin(surv_time, dropout_time)
  event         <- as.integer(surv_time <= dropout_time)
  calendar_time <- accrual_time + tte

  # Expand cell index into one column per factor
  sub_cols <- lapply(seq_along(sub_names),
                     function(f) as.integer(level_table[cell, f]))
  names(sub_cols) <- sub_names

  base_cols <- c(
    list(sim = rep(seq_len(nsim), each = n), group = group_id),
    sub_cols,
    list(accrual_time  = accrual_time,
         surv_time     = surv_time,
         dropout_time  = dropout_time,
         tte           = tte,
         event         = event,
         calendar_time = calendar_time)
  )

  data.frame(base_cols, stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------ #
#  Internal helper: piecewise exponential sampler
# ------------------------------------------------------------------ #
rpiece_exp_r <- function(n, hazard, e.time) {
  if (length(hazard) == 1L) {
    return(dqrng::dqrexp(n, rate = hazard))
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

  # C++ piecewise exponential sampler
  rpiece_exp_cpp(n, hazard, fin_time, cum_haz)
}

# ------------------------------------------------------------------ #
#  Internal helper: resolve e.time / d.time for a specific group
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
