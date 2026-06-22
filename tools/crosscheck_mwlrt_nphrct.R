# Cross-check survdiff_fast(weight = "mwlrt") against nphRCT::wlrt(method = "mw").
#
# Purpose: confirm the corrected modestly-weighted weight cap, 1 / S(t_star)
# (the pooled Kaplan-Meier value at t_star), matches the canonical nphRCT
# implementation in settings where the cap binds, i.e. t_star inside the
# observed event range so that events occur after t_star. One t_star is placed
# strictly between two event times to probe the boundary convention for
# S(t_star); on-event and between-event t_star giving the same agreement
# confirms the right-continuous pooled Kaplan-Meier convention is shared.
#
# Run with the package loaded (devtools::load_all() or library(FastSurvival)).
# Results are printed to the console; the script writes no files.

library(FastSurvival)

stopifnot(requireNamespace("nphRCT", quietly = TRUE),
          requireNamespace("survival", quietly = TRUE))

set.seed(20240622)

# ---- Two-group data set with a delayed separation ------------------------
# Group 1 = control, group 2 = experimental (simdata_fast convention).
# Continuous times keep event times distinct so t_star can sit between events.
sim <- simdata_fast(
  nsim = 1, n = c(250, 250),
  a.time = c(0, 6), a.prop = 1,
  e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 24)),
  e.time = list(NULL, c(0, 6, Inf)),
  d.median = list(60, 60),
  seed = 20240622
)

time  <- sim$tte
event <- sim$event
grp   <- sim$group
ctrl  <- 1L
arm01 <- as.integer(grp != ctrl)   # 1 = experimental

# ---- t_star grid: on an event time, between events, and a later point ----
et    <- sort(unique(time[event == 1L]))
i_mid <- as.integer(round(length(et) * 0.40))
t_on      <- et[i_mid]                            # exactly an event time
t_between <- (et[i_mid] + et[i_mid + 1L]) / 2     # strictly between two events
t_later   <- et[as.integer(round(length(et) * 0.65))]
t_grid    <- c(on_event = t_on, between_events = t_between, later = t_later)

# ---- nphRCT standardized statistic z = u / sqrt(v_u) ---------------------
get_num <- function(res, nm) {
  if (!is.null(res[[nm]])) as.numeric(res[[nm]]) else NA_real_
}

nphrct_z <- function(time, event, arm01, ts, show_str = FALSE) {
  df  <- data.frame(time = time, status = event, arm = arm01)
  res <- nphRCT::wlrt(survival::Surv(time, status) ~ arm, data = df,
                      method = "mw", t_star = ts)
  if (show_str) {
    cat("\n--- str() of one nphRCT::wlrt result (for field verification) ---\n")
    utils::str(res)
    cat("-----------------------------------------------------------------\n")
  }
  u  <- get_num(res, "u")
  vu <- get_num(res, "v_u")
  if (is.na(vu)) vu <- get_num(res, "v")   # fall back if the field is named 'v'
  u / sqrt(vu)
}

# ---- Compare -------------------------------------------------------------
rows <- lapply(seq_along(t_grid), function(k) {
  ts   <- t_grid[[k]]
  z_fs <- as.numeric(survdiff_fast(time, event, grp, control = ctrl,
                                   side = 1, weight = "mwlrt", t_star = ts))
  z_np <- nphrct_z(time, event, arm01, ts, show_str = (k == 1L))
  data.frame(
    case            = names(t_grid)[k],
    t_star          = ts,
    z_survdiff_fast = z_fs,
    z_nphRCT        = z_np,
    abs_diff        = abs(abs(z_fs) - abs(z_np))
  )
})
out <- do.call(rbind, rows)
out$match_1e_minus_6 <- out$abs_diff < 1e-6

cat("\n=== survdiff_fast (mwlrt) vs nphRCT::wlrt (mw) ===\n")
print(out, row.names = FALSE, digits = 10)

cat("\nNote: signs may differ by group-coding convention; the cap fix is",
    "\nconfirmed by agreement of |z| (abs_diff near 0). Agreement in the",
    "\n'between_events' case confirms the shared right-continuous S(t_star).\n",
    sep = "")

cat("\n")
print(sessionInfo())
