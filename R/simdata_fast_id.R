#' Illness-death (two-endpoint, switching) simulation
#'
#' Validates the transition-hazard arguments, builds the per-group
#' piecewise-exponential pre-computations, and calls simdata_core_id. The seed
#' is already set by the caller. With 'prevalence' the work is passed to
#' simdata_fast_id_sub(); without it the code below is unchanged.
#'
#' @inheritParams simdata_fast
#' @param alloc_given \code{TRUE} when \code{alloc} was given in the call of
#'   \code{simdata_fast()}.
#' @return The data frame returned by \code{simdata_fast()} for the
#'   illness-death model.
#' @noRd
simdata_fast_id <- function(nsim, n, alloc, a.time, a.rate, a.prop,
                            h01.hazard, h01.median, h01.time,
                            h02.hazard, h02.median, h02.time,
                            h12.hazard, h12.median, h12.time,
                            switch.prop,
                            h12.switch.hazard, h12.switch.median, h12.switch.time,
                            switch.clock,
                            d.hazard, d.median, d.time,
                            prevalence, fixed.alloc = FALSE,
                            alloc_given = FALSE,
                            h1s.hazard = NULL, h1s.median = NULL,
                            h1s.time = NULL,
                            h12.switch.clock = "switch",
                            switch.aft.factor = NULL,
                            h0s.hazard = NULL, h0s.median = NULL,
                            h0s.time = NULL,
                            h01.switch.hazard = NULL,
                            h01.switch.median = NULL,
                            h01.switch.time = NULL,
                            h02.switch.hazard = NULL,
                            h02.switch.median = NULL,
                            h02.switch.time = NULL,
                            e1.censor.switch = FALSE) {

  if (!(identical(switch.clock, "reset") ||
        identical(switch.clock, "forward"))) {
    stop("'switch.clock' must be \"reset\" or \"forward\"")
  }
  if (!(identical(h12.switch.clock, "switch") ||
        identical(h12.switch.clock, "same"))) {
    stop("'h12.switch.clock' must be \"switch\" or \"same\"")
  }
  if (!is.logical(e1.censor.switch) || length(e1.censor.switch) != 1L ||
      is.na(e1.censor.switch)) {
    stop("'e1.censor.switch' must be TRUE or FALSE")
  }

  # Resolve transition hazards (exactly one of hazard or median per quantity).
  if (!is.null(h01.hazard) && !is.null(h01.median)) {
    stop("Specify exactly one of 'h01.hazard' and 'h01.median'")
  }
  if (is.null(h01.hazard) && is.null(h01.median)) {
    stop("One of 'h01.hazard' or 'h01.median' must be supplied for the ",
         "illness-death model")
  }
  if (!is.null(h01.median)) h01.hazard <- convert_median_to_hazard(h01.median)

  if (!is.null(h02.hazard) && !is.null(h02.median)) {
    stop("Specify exactly one of 'h02.hazard' and 'h02.median'")
  }
  if (is.null(h02.hazard) && is.null(h02.median)) {
    stop("One of 'h02.hazard' or 'h02.median' must be supplied for the ",
         "illness-death model")
  }
  if (!is.null(h02.median)) h02.hazard <- convert_median_to_hazard(h02.median)

  if (!is.null(h12.hazard) && !is.null(h12.median)) {
    stop("Specify exactly one of 'h12.hazard' and 'h12.median'")
  }
  if (!is.null(h12.median)) h12.hazard <- convert_median_to_hazard(h12.median)

  if (!is.null(h12.switch.hazard) && !is.null(h12.switch.median)) {
    stop("Specify exactly one of 'h12.switch.hazard' and 'h12.switch.median'")
  }
  if (!is.null(h12.switch.median)) {
    h12.switch.hazard <- convert_median_to_hazard(h12.switch.median)
  }

  # Hazards of the switching rules (at most one of hazard and median each)
  one_of <- function(hz, md, nm) {
    if (!is.null(hz) && !is.null(md)) {
      stop("Specify exactly one of '", nm, ".hazard' and '", nm, ".median'")
    }
    if (!is.null(md)) convert_median_to_hazard(md) else hz
  }
  h1s.hazard <- one_of(h1s.hazard, h1s.median, "h1s")
  h0s.hazard <- one_of(h0s.hazard, h0s.median, "h0s")
  h01.switch.hazard <- one_of(h01.switch.hazard, h01.switch.median,
                              "h01.switch")
  h02.switch.hazard <- one_of(h02.switch.hazard, h02.switch.median,
                              "h02.switch")
  if (is.null(h0s.hazard) &&
      (!is.null(h01.switch.hazard) || !is.null(h02.switch.hazard))) {
    stop("'h01.switch.*' and 'h02.switch.*' are used only with 'h0s.*' ",
         "(switching before the intermediate event)")
  }
  if (isTRUE(e1.censor.switch) && is.null(h0s.hazard)) {
    stop("'e1.censor.switch' applies to switching before the intermediate ",
         "event, which needs 'h0s.*'")
  }
  if (!is.null(switch.aft.factor) && !is.null(h12.switch.hazard)) {
    stop("Specify at most one of 'switch.aft.factor' and ",
         "'h12.switch.hazard' (or 'h12.switch.median')")
  }
  new_rules <- !is.null(h1s.hazard) || !is.null(h0s.hazard) ||
    !is.null(switch.aft.factor) || identical(switch.clock, "forward") ||
    identical(h12.switch.clock, "same") || isTRUE(e1.censor.switch)

  has_dropout <- !is.null(d.hazard) || !is.null(d.median)
  if (!is.null(d.hazard) && !is.null(d.median)) {
    stop("Specify exactly one of 'd.hazard' and 'd.median'")
  }
  if (!is.null(d.median)) d.hazard <- convert_median_to_hazard(d.median)

  # The switching rules added in 2.0.0 are generated by the subgroup kernel.
  # Without 'prevalence' each group is passed to it as one cell, which draws
  # the same random numbers as simdata_core_id(), and the subgroup column is
  # dropped.
  sub_prev  <- prevalence
  sub_n     <- n
  sub_alloc <- alloc_given
  if (is.null(prevalence) && new_rules) {
    two <- (length(n) == 2L) || alloc_given ||
      is.list(h01.hazard) || is.list(h02.hazard) || is.list(h12.hazard) ||
      is.list(h12.switch.hazard) || is.list(switch.prop) ||
      is.list(d.hazard) || is.list(h1s.hazard) || is.list(h0s.hazard) ||
      is.list(h01.switch.hazard) || is.list(h02.switch.hazard) ||
      is.list(switch.aft.factor)
    if (length(n) == 1L && two) sub_n <- split_total(n, alloc)
    sub_prev  <- 1
    sub_alloc <- FALSE
  }
  if (!is.null(sub_prev)) {
    out <- simdata_fast_id_sub(
      nsim = nsim, n = sub_n, alloc = alloc, alloc_given = sub_alloc,
      a.time = a.time, a.rate = a.rate, a.prop = a.prop,
      h01.hazard = h01.hazard, h01.time = h01.time,
      h02.hazard = h02.hazard, h02.time = h02.time,
      h12.hazard = h12.hazard, h12.time = h12.time,
      switch.prop = switch.prop,
      h12.switch.hazard = h12.switch.hazard,
      h12.switch.time = h12.switch.time,
      has_dropout = has_dropout, d.hazard = d.hazard, d.time = d.time,
      prevalence = sub_prev, fixed.alloc = fixed.alloc,
      h1s.hazard = h1s.hazard, h1s.time = h1s.time,
      h0s.hazard = h0s.hazard, h0s.time = h0s.time,
      h01.switch.hazard = h01.switch.hazard,
      h01.switch.time = h01.switch.time,
      h02.switch.hazard = h02.switch.hazard,
      h02.switch.time = h02.switch.time,
      switch.aft.factor = switch.aft.factor, switch.clock = switch.clock,
      h12.switch.clock = h12.switch.clock,
      e1.censor.switch = e1.censor.switch)
    if (is.null(prevalence)) out <- out[names(out) != "subgroup"]
    return(out)
  }

  # A two-group simulation is signalled by a length-two 'n', an explicitly
  # supplied 'alloc', or any group-specific (length-two list) hazard, switch, or
  # dropout argument.
  is_two_group <- (length(n) == 2L) || alloc_given ||
    is.list(h01.hazard) || is.list(h02.hazard) ||
    is.list(h12.hazard) || is.list(h12.switch.hazard) ||
    is.list(switch.prop) || is.list(d.hazard)

  if (length(n) == 1L) {
    n_grp <- if (is_two_group) split_total(n, alloc) else n
  } else if (length(n) == 2L) {
    n_grp <- n
  } else {
    stop("'n' must be a scalar (total N) or a vector of length 2 (per-group)")
  }
  n_groups  <- if (is_two_group) 2L else 1L
  n_grp_int <- if (n_groups == 1L) as.integer(n_grp[1L]) else as.integer(n_grp[1:2])
  check_output_size(nsim, n_grp_int)
  # As for 'n', a whole number up to rounding error (such as 100 * 0.29) is
  # rounded rather than truncated.
  nsim <- round(nsim)
  g2        <- (n_groups == 2L)

  # As with subgroups (simdata_fast_id_sub()), a list is per group and must
  # have one element per group.
  if (g2) {
    spec_args <- list(h01.hazard = h01.hazard, h02.hazard = h02.hazard,
                      h12.hazard = h12.hazard,
                      h12.switch.hazard = h12.switch.hazard,
                      switch.prop = switch.prop, d.hazard = d.hazard)
    for (nm in names(spec_args)) {
      x <- spec_args[[nm]]
      if (is.list(x) && length(x) != 2L) {
        stop("For a two-group simulation, ", spec_arg_label(nm),
             " must be a single ",
             "(shared) specification or a list of length 2 (one element per ",
             "group); it is a list of length ", length(x), ".")
      }
    }
  }

  # Resolve a possibly group-specific argument for group g (1 = control,
  # 2 = treatment). A length-two list is per group; anything else is shared.
  id_grp <- function(x, g) if (is.list(x)) x[[g]] else x

  get_pi <- function(g) {
    if (is.null(switch.prop)) return(0)
    p <- id_grp(switch.prop, g)
    if (is.null(p)) return(0)
    if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 0 ||
        p > 1) {
      stop("Each 'switch.prop' value must be a single probability in [0, 1]")
    }
    as.numeric(p)
  }
  pi_c <- get_pi(1L)
  pi_t <- if (g2) get_pi(2L) else 0

  if ((pi_c > 0 || pi_t > 0) && is.null(h12.switch.hazard)) {
    stop("'h12.switch.hazard' (or 'h12.switch.median') must be supplied when ",
         "any 'switch.prop' is positive")
  }

  # Post-event no-switch hazard defaults to the direct terminal hazard h02.
  get_h12      <- function(g) if (is.null(h12.hazard)) id_grp(h02.hazard, g) else id_grp(h12.hazard, g)
  get_h12_time <- function(g) if (is.null(h12.hazard)) id_grp(h02.time, g)   else id_grp(h12.time, g)

  h12_nm <- if (is.null(h12.hazard)) "h02" else "h12"
  h01_c <- piecewise_precompute(id_grp(h01.hazard, 1L), id_grp(h01.time, 1L), "h01")
  h01_t <- if (g2) piecewise_precompute(id_grp(h01.hazard, 2L), id_grp(h01.time, 2L), "h01") else h01_c
  h02_c <- piecewise_precompute(id_grp(h02.hazard, 1L), id_grp(h02.time, 1L), "h02")
  h02_t <- if (g2) piecewise_precompute(id_grp(h02.hazard, 2L), id_grp(h02.time, 2L), "h02") else h02_c
  h12_c <- piecewise_precompute(get_h12(1L), get_h12_time(1L), h12_nm)
  h12_t <- if (g2) piecewise_precompute(get_h12(2L), get_h12_time(2L), h12_nm) else h12_c

  empty_spec <- list(hazard = 1, fin_time = numeric(0), cum_haz = numeric(0))
  if (pi_c > 0 || pi_t > 0) {
    h12s_c <- piecewise_precompute(id_grp(h12.switch.hazard, 1L), id_grp(h12.switch.time, 1L), "h12.switch")
    h12s_t <- if (g2) {
      piecewise_precompute(id_grp(h12.switch.hazard, 2L), id_grp(h12.switch.time, 2L), "h12.switch")
    } else {
      h12s_c
    }
  } else {
    h12s_c <- empty_spec
    h12s_t <- empty_spec
  }

  if (has_dropout) {
    d_c <- piecewise_precompute(id_grp(d.hazard, 1L), id_grp(d.time, 1L), "d")
    d_t <- if (g2) piecewise_precompute(id_grp(d.hazard, 2L), id_grp(d.time, 2L), "d") else d_c
  } else {
    d_c <- empty_spec
    d_t <- empty_spec
  }

  acc <- resolve_accrual_counts(n_grp_int, a.time, a.rate, a.prop)

  simdata_core_id(
    as.integer(nsim), n_grp_int,
    as.numeric(acc$a.time_full), acc$acc_counts_c, acc$acc_counts_t,
    h01_c$hazard, h01_c$fin_time, h01_c$cum_haz,
    h01_t$hazard, h01_t$fin_time, h01_t$cum_haz,
    h02_c$hazard, h02_c$fin_time, h02_c$cum_haz,
    h02_t$hazard, h02_t$fin_time, h02_t$cum_haz,
    h12_c$hazard, h12_c$fin_time, h12_c$cum_haz,
    h12_t$hazard, h12_t$fin_time, h12_t$cum_haz,
    h12s_c$hazard, h12s_c$fin_time, h12s_c$cum_haz,
    h12s_t$hazard, h12s_t$fin_time, h12s_t$cum_haz,
    as.numeric(pi_c), as.numeric(pi_t),
    has_dropout,
    d_c$hazard, d_c$fin_time, d_c$cum_haz,
    d_t$hazard, d_t$fin_time, d_t$cum_haz
  )
}
