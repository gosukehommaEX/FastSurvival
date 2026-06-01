#' Fast Sequential Analysis of Simulated Trial Data
#'
#' @description
#' Performs interim or sequential analyses of simulated two-group
#' time-to-event data at one or more analysis times ("looks"). Each look is
#' defined either by a target cumulative number of events
#' (information-based timing) or by a calendar time (calendar-based timing).
#' At every look the data are administratively censored at the corresponding
#' calendar cutoff, and the requested statistics are computed for each
#' simulated trial by reusing \code{\link{survdiff_fast}},
#' \code{\link{coxph_fast}}, \code{\link{rmst_fast}}, \code{\link{survfit_fast}},
#' \code{\link{maxcombo_fast}}, and \code{\link{ahsw_fast}}. Optionally the same
#' statistics are also reported within each subgroup. The censoring and time
#' sorting are handled by a C++ backend, and the analysis cores are called with
#' \code{presorted = TRUE}, so each look avoids a redundant sort.
#'
#' @details
#' The input \code{data} is the data frame returned by
#' \code{\link{simdata_fast}} for a two-group trial. The columns \code{sim},
#' \code{group}, \code{accrual_time}, \code{tte}, and \code{event} are
#' required.
#'
#' For a look at calendar time \code{cutoff}, each subject with accrual time
#' \code{a} contributes only if enrolled by then (\code{a <= cutoff}). The
#' observed time at the look is \code{min(tte, cutoff - a)} and the observed
#' event indicator is the original \code{event} when the natural event or
#' dropout occurred on or before \code{cutoff} (\code{a + tte <= cutoff}), and
#' zero otherwise (administrative censoring at the look).
#'
#' When \code{event.looks} is supplied, the calendar cutoff for a target of
#' \code{d} events is the calendar time of the \code{d}-th event in that
#' simulated trial, counted over the whole trial population. If a simulation
#' contains fewer than \code{d} events, the target is never reached: the full
#' data are used, \code{reached} is \code{FALSE}, and \code{cutoff} is
#' \code{NA}. When \code{time.looks} is supplied, the cutoff is the specified
#' calendar time and \code{reached} is always \code{TRUE}. In both cases the
#' cutoff is determined once on the whole population and then used for the
#' overall analysis and for every subgroup analysis at that look.
#'
#' Exactly one of \code{event.looks} and \code{time.looks} must be supplied.
#'
#' The statistics are selected with \code{stat}, which may name one or more of
#' \code{"logrank"}, \code{"coxph"}, \code{"rmst"}, \code{"km"},
#' \code{"maxcombo"}, and \code{"ahsw"}.
#'
#' The \code{"logrank"} statistic is configurable. By default it is the
#' ordinary unweighted, unstratified two-group log-rank test and reproduces the
#' behaviour of earlier versions of this function exactly. A non-default
#' \code{weight} selects a weighted log-rank test (Fleming-Harrington,
#' modestly-weighted, Gehan-Breslow, or Tarone-Ware) for non-proportional
#' hazards, and a non-\code{NULL} \code{strata} selects the stratified test,
#' summing the per-stratum contributions. The two options combine to give the
#' stratified weighted log-rank test. Whatever configuration is chosen, the
#' result is written to the same \code{logrank.z}, \code{logrank.chisq}, and
#' \code{logrank.p} columns; the unweighted unstratified case is the \code{K = 1}
#' single-stratum degenerate form of the general statistic. The columns
#' \code{rho}, \code{gamma}, and \code{t_star} parametrize \code{weight}, and
#' \code{strata} names one or more subgroup columns of \code{data} used as the
#' stratification variable. The stratification is determined on the whole cut
#' data, independently of the \code{population} marginalization, so a stratified
#' overall analysis is the canonical primary test; in a single-subgroup
#' population the stratum is constant and the stratified test degenerates to the
#' ordinary one within that subset.
#'
#' The \code{"maxcombo"} statistic is the max-combo test of \code{mc.rho} and
#' \code{mc.gamma} Fleming-Harrington weights. Its \code{maxcombo.stat} is the
#' most extreme component (\code{min} of the component Z-scores when
#' \code{side = 1}, so that a negative value favours treatment, and the maximum
#' absolute component when \code{side = 2}), and \code{maxcombo.p} is the joint
#' multivariate-normal p-value, which already follows \code{side}.
#'
#' The \code{"ahsw"} statistic is the average hazard with survival weight of Uno
#' and Horiguchi on the window from 0 to \code{tau}. It reports the per-group
#' average hazards, the ratio (RAH) and difference (DAH) contrasts with their
#' confidence intervals, and two-sided p-values for both contrasts. The AHSW
#' p-values are always two-sided and do not depend on \code{side}, matching
#' \code{\link{ahsw_fast}}.
#'
#' When \code{by.subgroup = TRUE}, the output is given in long form with a
#' \code{population} column. Each \code{(sim, look)} produces one row for the
#' whole trial (\code{population = "overall"}) plus one row per subgroup level
#' of each subgroup factor in \code{data}. Subgroup factors are the columns
#' named \code{subgroup} or \code{subgroup1}, \code{subgroup2}, and so on. Each
#' factor is marginalized separately, so a factor with levels 1 and 2 yields
#' the populations \code{subgroup_1} and \code{subgroup_2} (or
#' \code{subgroup1_1}, \code{subgroup1_2}, and so on for numbered factors).
#' The look cutoff is always determined on the whole population, so subgroup
#' rows at a given look share the same \code{cutoff}, \code{reached}, and
#' \code{look.value}; their \code{n.enrolled} and \code{n.event} are the counts
#' within that subgroup. When \code{by.subgroup = FALSE} (default) the output
#' has no \code{population} column and one row per \code{(sim, look)}, matching
#' the whole-population analysis.
#'
#' @param data A data frame from \code{\link{simdata_fast}} for a two-group
#'   trial, containing at least \code{sim}, \code{group}, \code{accrual_time},
#'   \code{tte}, and \code{event}.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group.
#' @param event.looks A numeric vector of target cumulative event counts, one
#'   per look. Mutually exclusive with \code{time.looks}.
#' @param time.looks A numeric vector of calendar times, one per look.
#'   Mutually exclusive with \code{event.looks}.
#' @param stat A character vector naming the statistics to compute. Any subset
#'   of \code{"logrank"}, \code{"coxph"}, \code{"rmst"}, \code{"km"},
#'   \code{"maxcombo"}, and \code{"ahsw"}. Defaults to \code{"logrank"}.
#' @param tau A single positive numeric value, the restriction horizon for
#'   \code{"rmst"} and the truncation time for \code{"ahsw"}. Required only when
#'   \code{"rmst"} or \code{"ahsw"} is requested.
#' @param t.eval A single positive numeric value, the landmark time for
#'   \code{"km"}. Required only when \code{"km"} is requested.
#' @param conf.int A single numeric value in (0, 1), the confidence level for
#'   \code{"coxph"}, \code{"rmst"}, and \code{"ahsw"}. Defaults to 0.95.
#' @param side An integer, either 1 or 2. When \code{side = 2} (default),
#'   two-sided p-values \code{2 pnorm(-|z|)} are reported for log-rank, Cox, and
#'   RMST, and the two-sided max-combo test is used. When \code{side = 1}, the
#'   one-sided p-value in the direction of treatment benefit is reported for
#'   each of those statistics, and the one-sided max-combo test is used. The
#'   test statistics themselves are always reported with their natural sign, so
#'   the choice of \code{side} affects only the p-value columns. For log-rank
#'   and Cox the benefit direction is a negative Z (the one-sided p-value is the
#'   lower tail \code{pnorm(z)}); for RMST it is a positive Z (the upper tail
#'   \code{pnorm(-z)}). The AHSW p-values are always two-sided. For
#'   group-sequential boundary comparisons (for example with gsDesign or
#'   rpact), align the sign of the reported Z with the boundary convention
#'   before comparing.
#' @param by.subgroup A logical value. When \code{TRUE}, the analysis is also
#'   reported within each subgroup level, and the output gains a
#'   \code{population} column (long form). When \code{FALSE} (default), only
#'   the whole-population analysis is returned and no \code{population} column
#'   is added. Requires at least one subgroup column in \code{data}.
#' @param weight A character string naming the weight scheme for the
#'   \code{"logrank"} statistic. \code{"logrank"} (default) is the ordinary
#'   unweighted test. \code{"fh"} is the Fleming-Harrington G(rho, gamma) test,
#'   \code{"mwlrt"} is the modestly-weighted log-rank test (requires
#'   \code{t_star}), \code{"gehan"} is Gehan-Breslow, and \code{"tarone-ware"}
#'   is Tarone-Ware. See \code{\link{survdiff_fast}}.
#' @param rho A numeric Fleming-Harrington first parameter, used only when
#'   \code{weight = "fh"}. Defaults to 0.
#' @param gamma A numeric Fleming-Harrington second parameter, used only when
#'   \code{weight = "fh"}. Defaults to 0.
#' @param t_star A single non-negative numeric value, the timepoint of the
#'   modestly-weighted log-rank test. Required only when \code{weight =
#'   "mwlrt"}.
#' @param strata An optional character vector naming one or more subgroup
#'   columns of \code{data} to use as the stratification variable for the
#'   \code{"logrank"} statistic. \code{NULL} (default) gives the unstratified
#'   test. When several columns are named their interaction defines the strata.
#'   Stratification applies only to the \code{"logrank"} statistic; the other
#'   statistics ignore it.
#' @param mc.rho A numeric vector of Fleming-Harrington first parameters for the
#'   \code{"maxcombo"} statistic, one per component weight. Defaults to
#'   \code{c(0, 0, 1, 1)}.
#' @param mc.gamma A numeric vector of Fleming-Harrington second parameters for
#'   the \code{"maxcombo"} statistic, aligned with \code{mc.rho}. Defaults to
#'   \code{c(0, 1, 0, 1)}.
#' @param abseps A single positive numeric value, the absolute error tolerance
#'   passed to the multivariate normal integration of the \code{"maxcombo"}
#'   p-value. Defaults to 1e-5.
#' @param maxpts A single positive integer, the maximum number of function
#'   evaluations for the quasi-Monte-Carlo integration used by the
#'   \code{"maxcombo"} p-value when four or more weights are supplied. Defaults
#'   to 25000.
#'
#' @return A data frame. When \code{by.subgroup = FALSE}, it has
#'   \code{nsim * length(looks)} rows. When \code{by.subgroup = TRUE}, it has
#'   \code{nsim * length(looks) * (1 + total subgroup levels)} rows and an
#'   extra \code{population} column placed after \code{look.value}. The common
#'   columns are \code{sim}, \code{look} (1-based look index), \code{look.type}
#'   (\code{"event"} or \code{"time"}), \code{look.value} (the requested event
#'   count or calendar time), optionally \code{population}, \code{cutoff} (the
#'   calendar time used, \code{NA} when an event target was not reached),
#'   \code{reached}, \code{n.enrolled}, and \code{n.event}, followed by the
#'   columns of the requested statistics. A statistic that cannot be computed
#'   for a row (no events, or an empty group) is \code{NA}. The statistic
#'   columns are \code{logrank.z}, \code{logrank.chisq}, and \code{logrank.p}
#'   for \code{"logrank"}; \code{cox.coef}, \code{cox.hr}, \code{cox.se},
#'   \code{cox.z}, \code{cox.p}, \code{cox.lower}, and \code{cox.upper} for
#'   \code{"coxph"}; \code{rmst.ctrl}, \code{rmst.trt}, \code{rmst.diff},
#'   \code{rmst.diff.lower}, \code{rmst.diff.upper}, \code{rmst.z}, and
#'   \code{rmst.p} for \code{"rmst"}; \code{km.surv.ctrl} and
#'   \code{km.surv.trt} for \code{"km"}; \code{maxcombo.stat} and
#'   \code{maxcombo.p} for \code{"maxcombo"}; and \code{ahsw.ah.ctrl},
#'   \code{ahsw.ah.trt}, \code{ahsw.rah}, \code{ahsw.rah.lower},
#'   \code{ahsw.rah.upper}, \code{ahsw.p.rah}, \code{ahsw.dah},
#'   \code{ahsw.dah.lower}, \code{ahsw.dah.upper}, and \code{ahsw.p.dah} for
#'   \code{"ahsw"}. The Z columns \code{logrank.z}, \code{cox.z}, and
#'   \code{rmst.z} carry the natural sign of each test, and the p-value columns
#'   follow \code{side} except for the AHSW p-values, which are two-sided.
#'
#' @examples
#' df <- simdata_fast(
#'   nsim       = 50,
#'   n          = c(150, 150),
#'   a.time     = c(0, 12),
#'   a.rate     = 300 / 12,
#'   e.hazard   = list(list(0.10, 0.07), 0.05),
#'   prevalence = c(0.5, 0.5),
#'   seed       = 1
#' )
#'
#' # Whole-population analysis at two event-based looks
#' res1 <- analysis_fast(df, control = 1, event.looks = c(80, 140))
#' head(res1)
#'
#' # Stratified log-rank on the subgroup factor
#' res2 <- analysis_fast(df, control = 1, time.looks = 24,
#'                       stat = "logrank", strata = "subgroup")
#' head(res2)
#'
#' # Fleming-Harrington G(0, 1) weighted log-rank for delayed effects
#' res3 <- analysis_fast(df, control = 1, time.looks = 24,
#'                       stat = "logrank", weight = "fh", rho = 0, gamma = 1)
#' head(res3)
#'
#' # Max-combo and AHSW together
#' res4 <- analysis_fast(df, control = 1, time.looks = 24,
#'                       stat = c("maxcombo", "ahsw"), tau = 18)
#' head(res4)
#'
#' @seealso \code{\link{simdata_fast}}, \code{\link{survdiff_fast}},
#'   \code{\link{coxph_fast}}, \code{\link{rmst_fast}},
#'   \code{\link{survfit_fast}}, \code{\link{maxcombo_fast}},
#'   \code{\link{ahsw_fast}}.
#'
#' @importFrom stats pnorm
#' @export
analysis_fast <- function(data, control,
                          event.looks = NULL, time.looks = NULL,
                          stat = "logrank",
                          tau = NULL, t.eval = NULL,
                          conf.int = 0.95,
                          side = 2,
                          by.subgroup = FALSE,
                          weight = c("logrank", "fh", "mwlrt", "gehan",
                                     "tarone-ware"),
                          rho = 0, gamma = 0, t_star = NULL,
                          strata = NULL,
                          mc.rho = c(0, 0, 1, 1), mc.gamma = c(0, 1, 0, 1),
                          abseps = 1e-5, maxpts = 25000) {

  weight <- match.arg(weight)

  # ------------------------------------------------------------------ #
  #  Validate inputs
  # ------------------------------------------------------------------ #
  req_cols <- c("sim", "group", "accrual_time", "tte", "event")
  if (!is.data.frame(data) || !all(req_cols %in% names(data))) {
    stop("'data' must be a data frame with columns: ",
         paste(req_cols, collapse = ", "))
  }
  if (missing(control)) {
    stop("'control' must be supplied")
  }

  has_event <- !is.null(event.looks)
  has_time  <- !is.null(time.looks)
  if (has_event == has_time) {
    stop("Supply exactly one of 'event.looks' and 'time.looks'")
  }
  look.type   <- if (has_event) "event" else "time"
  look.values <- if (has_event) event.looks else time.looks
  if (length(look.values) < 1L || any(!is.finite(look.values)) ||
      any(look.values <= 0)) {
    stop("Look values must be positive and finite")
  }
  if (has_event && any(look.values != round(look.values))) {
    stop("'event.looks' must contain whole numbers (event counts)")
  }
  n.looks <- length(look.values)

  allowed_stat <- c("logrank", "coxph", "rmst", "km", "maxcombo", "ahsw")
  if (!all(stat %in% allowed_stat)) {
    stop("'stat' must be a subset of: ",
         paste(allowed_stat, collapse = ", "))
  }
  stat <- unique(stat)

  if (conf.int <= 0 || conf.int >= 1) {
    stop("'conf.int' must be in (0, 1)")
  }
  if (length(side) != 1L || !side %in% c(1, 2)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided)")
  }
  side <- as.integer(side)
  if (("rmst" %in% stat || "ahsw" %in% stat) &&
      (is.null(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0)) {
    stop("'tau' must be a single positive value when 'rmst' or 'ahsw' is ",
         "requested")
  }
  if ("km" %in% stat &&
      (is.null(t.eval) || length(t.eval) != 1L || !is.finite(t.eval) ||
       t.eval <= 0)) {
    stop("'t.eval' must be a single positive value when 'km' is requested")
  }
  if (weight == "mwlrt" &&
      (is.null(t_star) || length(t_star) != 1L || !is.finite(t_star) ||
       t_star < 0)) {
    stop("'t_star' must be a single non-negative value when weight = 'mwlrt'")
  }
  if ("maxcombo" %in% stat && length(mc.rho) != length(mc.gamma)) {
    stop("'mc.rho' and 'mc.gamma' must have the same length")
  }

  # ------------------------------------------------------------------ #
  #  Subgroup columns and population definitions
  # ------------------------------------------------------------------ #
  sub_cols <- grep("^subgroup([0-9]+)?$", names(data), value = TRUE)
  if (length(sub_cols) > 1L) {
    sfx      <- suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", sub_cols)))
    sub_cols <- sub_cols[order(sfx, sub_cols)]
  }

  if (by.subgroup && length(sub_cols) == 0L) {
    stop("'by.subgroup = TRUE' but 'data' has no subgroup columns")
  }

  # Stratification columns (for the log-rank statistic). Independent of the
  # 'population' marginalization done by 'by.subgroup'.
  use_strata <- !is.null(strata)
  if (use_strata) {
    if (!is.character(strata) || length(strata) < 1L ||
        !all(strata %in% names(data))) {
      stop("'strata' must name one or more columns of 'data'")
    }
  }

  # Population list: 'overall' first, then one entry per (factor, level).
  # Each entry stores the column name, the level value, and the row label.
  pop_defs <- list(list(col = NA_character_, level = NA_integer_,
                        label = "overall"))
  if (by.subgroup) {
    for (cn in sub_cols) {
      levs <- sort(unique(data[[cn]]))
      for (lv in levs) {
        pop_defs[[length(pop_defs) + 1L]] <-
          list(col = cn, level = lv, label = paste0(cn, "_", lv))
      }
    }
  }
  n_pop <- length(pop_defs)

  # ------------------------------------------------------------------ #
  #  Treatment indicator, per-simulation row index, subgroup vectors
  # ------------------------------------------------------------------ #
  grp <- data$group
  if (is.factor(grp)) grp <- as.character(grp)
  j_all <- as.integer(grp != control)

  sg_all <- if (by.subgroup) lapply(sub_cols, function(cn) data[[cn]]) else NULL
  if (by.subgroup) names(sg_all) <- sub_cols

  # Stratum label per row: a single column as-is, or the interaction of several
  st_all <- if (use_strata) {
    if (length(strata) == 1L) {
      as.character(data[[strata]])
    } else {
      do.call(paste, c(lapply(strata, function(cn) data[[cn]]), sep = "."))
    }
  } else NULL

  idx_by_sim <- split(seq_len(nrow(data)), data$sim)
  ord_names  <- order(as.numeric(names(idx_by_sim)))
  idx_by_sim <- idx_by_sim[ord_names]
  sim_vals   <- as.numeric(names(idx_by_sim))
  nsim       <- length(idx_by_sim)
  total      <- nsim * n.looks * n_pop

  # ------------------------------------------------------------------ #
  #  Pre-allocate output columns
  # ------------------------------------------------------------------ #
  out <- list(
    sim        = numeric(total),
    look       = integer(total),
    look.type  = rep(look.type, total),
    look.value = numeric(total)
  )
  if (by.subgroup) out$population <- character(total)
  out$cutoff     <- numeric(total)
  out$reached    <- logical(total)
  out$n.enrolled <- integer(total)
  out$n.event    <- integer(total)

  if ("logrank" %in% stat) {
    out$logrank.z     <- rep(NA_real_, total)
    out$logrank.chisq <- rep(NA_real_, total)
    out$logrank.p     <- rep(NA_real_, total)
  }
  if ("coxph" %in% stat) {
    out$cox.coef  <- rep(NA_real_, total)
    out$cox.hr    <- rep(NA_real_, total)
    out$cox.se    <- rep(NA_real_, total)
    out$cox.z     <- rep(NA_real_, total)
    out$cox.p     <- rep(NA_real_, total)
    out$cox.lower <- rep(NA_real_, total)
    out$cox.upper <- rep(NA_real_, total)
  }
  if ("rmst" %in% stat) {
    out$rmst.ctrl       <- rep(NA_real_, total)
    out$rmst.trt        <- rep(NA_real_, total)
    out$rmst.diff       <- rep(NA_real_, total)
    out$rmst.diff.lower <- rep(NA_real_, total)
    out$rmst.diff.upper <- rep(NA_real_, total)
    out$rmst.z          <- rep(NA_real_, total)
    out$rmst.p          <- rep(NA_real_, total)
  }
  if ("km" %in% stat) {
    out$km.surv.ctrl <- rep(NA_real_, total)
    out$km.surv.trt  <- rep(NA_real_, total)
  }
  if ("maxcombo" %in% stat) {
    out$maxcombo.stat <- rep(NA_real_, total)
    out$maxcombo.p    <- rep(NA_real_, total)
  }
  if ("ahsw" %in% stat) {
    out$ahsw.ah.ctrl   <- rep(NA_real_, total)
    out$ahsw.ah.trt    <- rep(NA_real_, total)
    out$ahsw.rah       <- rep(NA_real_, total)
    out$ahsw.rah.lower <- rep(NA_real_, total)
    out$ahsw.rah.upper <- rep(NA_real_, total)
    out$ahsw.p.rah     <- rep(NA_real_, total)
    out$ahsw.dah       <- rep(NA_real_, total)
    out$ahsw.dah.lower <- rep(NA_real_, total)
    out$ahsw.dah.upper <- rep(NA_real_, total)
    out$ahsw.p.dah     <- rep(NA_real_, total)
  }

  # ------------------------------------------------------------------ #
  #  Two-sided / one-sided p-value from a signed Z statistic.
  #  'eff_dir' is the sign of Z that corresponds to treatment benefit:
  #  -1 for log-rank and Cox (benefit gives a negative Z), +1 for RMST
  #  (benefit gives a positive Z). The one-sided p-value is the tail in the
  #  benefit direction, so the three statistics share the same one-sided
  #  interpretation despite their differing natural signs.
  # ------------------------------------------------------------------ #
  p_from_z <- function(z, eff_dir) {
    if (is.na(z)) return(NA_real_)
    if (side == 2L) {
      2 * stats::pnorm(-abs(z))
    } else if (eff_dir < 0) {
      stats::pnorm(z)
    } else {
      stats::pnorm(-z)
    }
  }

  # ------------------------------------------------------------------ #
  #  Statistic writer for one population (subset of a look's cut data).
  #  'sc' is the stratum label vector aligned with (tc, ec, jc), or NULL when
  #  stratification is not requested.
  # ------------------------------------------------------------------ #
  write_stats <- function(pos, tc, ec, jc, sc) {
    n_ev <- sum(ec)
    n0   <- sum(jc == 0L)
    n1   <- sum(jc == 1L)
    both <- n0 > 0L && n1 > 0L

    out$n.enrolled[pos] <<- length(tc)
    out$n.event[pos]    <<- n_ev

    if ("logrank" %in% stat && n_ev > 0L && both) {
      # Unstratified: the cut data are already time-sorted, so presorted = TRUE.
      # Stratified: the cut data are sorted by time only, not by stratum then
      # time, so presorted = FALSE and survdiff_fast re-sorts internally.
      z <- if (is.null(sc)) {
        as.numeric(survdiff_fast(tc, ec, jc, control = 0L, side = 1L,
                                 presorted = TRUE, weight = weight,
                                 rho = rho, gamma = gamma, t_star = t_star))
      } else {
        as.numeric(survdiff_fast(tc, ec, jc, control = 0L, side = 1L,
                                 presorted = FALSE, strata = sc,
                                 weight = weight, rho = rho, gamma = gamma,
                                 t_star = t_star))
      }
      out$logrank.z[pos]     <<- z
      out$logrank.chisq[pos] <<- z * z
      out$logrank.p[pos]     <<- p_from_z(z, eff_dir = -1)
    }
    if ("coxph" %in% stat && n_ev > 0L && both) {
      cx <- coxph_fast(tc, ec, jc, control = 0L,
                       conf.int = conf.int, presorted = TRUE)
      coef_v <- unname(cx[1L])
      se_v   <- unname(cx[3L])
      z_v    <- if (is.finite(se_v) && se_v > 0) coef_v / se_v else NA_real_
      out$cox.coef[pos]  <<- coef_v
      out$cox.hr[pos]    <<- unname(cx[2L])
      out$cox.se[pos]    <<- se_v
      out$cox.z[pos]     <<- z_v
      out$cox.p[pos]     <<- p_from_z(z_v, eff_dir = -1)
      out$cox.lower[pos] <<- unname(cx[4L])
      out$cox.upper[pos] <<- unname(cx[5L])
    }
    if ("rmst" %in% stat && both) {
      rm <- rmst_fast(tc, ec, group = jc, control = 0L,
                      tau = tau, conf.int = conf.int, presorted = TRUE)
      z_d <- unname(rm["z.diff"])
      out$rmst.ctrl[pos]       <<- unname(rm["rmst.ctrl"])
      out$rmst.trt[pos]        <<- unname(rm["rmst.trt"])
      out$rmst.diff[pos]       <<- unname(rm["diff"])
      out$rmst.diff.lower[pos] <<- unname(rm["diff.lower"])
      out$rmst.diff.upper[pos] <<- unname(rm["diff.upper"])
      out$rmst.z[pos]          <<- z_d
      out$rmst.p[pos]          <<- p_from_z(z_d, eff_dir = 1)
    }
    if ("km" %in% stat) {
      if (n0 > 0L) {
        is0 <- jc == 0L
        kc  <- survfit_fast(tc[is0], ec[is0], t_eval = t.eval, presorted = TRUE)
        out$km.surv.ctrl[pos] <<- unname(kc["surv"])
      }
      if (n1 > 0L) {
        is1 <- jc == 1L
        kt  <- survfit_fast(tc[is1], ec[is1], t_eval = t.eval, presorted = TRUE)
        out$km.surv.trt[pos] <<- unname(kt["surv"])
      }
    }
    if ("maxcombo" %in% stat && n_ev > 0L && both) {
      mcf <- maxcombo_fast(tc, ec, jc, control = 0L, side = side,
                           rho = mc.rho, gamma = mc.gamma, presorted = TRUE,
                           abseps = abseps, maxpts = maxpts)
      out$maxcombo.stat[pos] <<- unname(mcf["statistic"])
      out$maxcombo.p[pos]    <<- unname(mcf["p.value"])
    }
    if ("ahsw" %in% stat && n_ev > 0L && both) {
      ah <- ahsw_fast(tc, ec, group = jc, control = 0L, tau = tau,
                      conf.int = conf.int, presorted = TRUE)
      out$ahsw.ah.ctrl[pos]   <<- unname(ah["ah.ctrl"])
      out$ahsw.ah.trt[pos]    <<- unname(ah["ah.trt"])
      out$ahsw.rah[pos]       <<- unname(ah["rah"])
      out$ahsw.rah.lower[pos] <<- unname(ah["rah.lower"])
      out$ahsw.rah.upper[pos] <<- unname(ah["rah.upper"])
      out$ahsw.p.rah[pos]     <<- unname(ah["p.rah"])
      out$ahsw.dah[pos]       <<- unname(ah["dah"])
      out$ahsw.dah.lower[pos] <<- unname(ah["dah.lower"])
      out$ahsw.dah.upper[pos] <<- unname(ah["dah.upper"])
      out$ahsw.p.dah[pos]     <<- unname(ah["p.dah"])
    }
  }

  # ------------------------------------------------------------------ #
  #  Main loop over simulations and looks
  # ------------------------------------------------------------------ #
  pos <- 0L
  for (si in seq_len(nsim)) {
    rows <- idx_by_sim[[si]]
    acc  <- data$accrual_time[rows]
    tt   <- data$tte[rows]
    ev   <- as.integer(data$event[rows])
    jj   <- j_all[rows]
    sg   <- if (by.subgroup) lapply(sg_all, function(v) v[rows]) else NULL
    st   <- if (use_strata) st_all[rows] else NULL

    # Calendar times of events (for event-driven cutoffs, whole population)
    cal_ev <- if (look.type == "event") (acc + tt)[ev == 1L] else NULL

    for (l in seq_len(n.looks)) {
      if (look.type == "event") {
        cutoff    <- nth_event_time_core(cal_ev, as.integer(look.values[l]))
        reached_l <- is.finite(cutoff)
      } else {
        cutoff    <- look.values[l]
        reached_l <- TRUE
      }

      cut    <- analysis_cut_core(acc, tt, ev, jj, cutoff)
      t_cut  <- cut$time
      e_cut  <- cut$event
      j_cut  <- cut$j
      idx    <- cut$idx
      cut_v  <- if (reached_l) cutoff else NA_real_

      # Subgroup labels and stratum labels realigned to the sorted cut order
      sg_cut <- if (by.subgroup) lapply(sg, function(v) v[idx]) else NULL
      st_cut <- if (use_strata) st[idx] else NULL

      for (pi in seq_len(n_pop)) {
        pos <- pos + 1L
        pd  <- pop_defs[[pi]]

        out$sim[pos]        <- sim_vals[si]
        out$look[pos]       <- l
        out$look.value[pos] <- look.values[l]
        out$cutoff[pos]     <- cut_v
        out$reached[pos]    <- reached_l
        if (by.subgroup) out$population[pos] <- pd$label

        if (is.na(pd$col)) {
          write_stats(pos, t_cut, e_cut, j_cut, st_cut)
        } else {
          keep <- sg_cut[[pd$col]] == pd$level
          sc_k <- if (use_strata) st_cut[keep] else NULL
          write_stats(pos, t_cut[keep], e_cut[keep], j_cut[keep], sc_k)
        }
      }
    }
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}
