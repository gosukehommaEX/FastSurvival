# Agreement and speed of the FastSurvival analysis functions against reference
# implementations on CRAN.
#
# Part 1 compares each function with its reference on the gbsg data of the
# survival package (686 patients, recurrence-free survival time in days,
# hormonal therapy as the group), as in the validation vignette, and records
# the largest absolute difference of the compared quantities. The window mean
# survival time and the weighted Kaplan-Meier test, which have no reference on
# CRAN, are compared with a Kaplan-Meier integral from survfit() and with the
# restricted mean survival time identity. The simulation layer is checked by
# comparing cutoff_fast() with simtrial::get_analysis_date(), and the log-rank
# and RMST statistics of analysis_fast() with survival::survdiff() and
# survRM2::rmst2() applied to the data censored at the same cutoffs.
#
# Part 2 times each function against its reference on a simulated two-group
# data set of 500 subjects with microbenchmark (1000 replicates), with the
# pre-sorted fast path and one-sided tests where applicable, as in
# tools/benchmark_speed.R. Only references on CRAN are timed.
#
# Results are written to tools/paper/output/bench_functions_agreement.csv,
# bench_functions_speed.csv, and bench_functions_data.csv.
#
# Run from the package root after installing FastSurvival 1.1.0 from CRAN
# (checked by machine_info.R):
#   source("tools/paper/data/bench_functions.R")

library(FastSurvival)
source(file.path("tools", "paper", "data", "machine_info.R"))

for (p in c("survival", "survRM2", "nph", "nphRCT", "survAH", "simtrial",
            "microbenchmark")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Package '", p, "' is required for this script.", call. = FALSE)
  }
}
library(survival)

out_dir <- file.path("tools", "paper", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Part 1: agreement with the references ----------------------------------
agree <- list()
add_agree <- function(method, reference, quantity, fast, ref) {
  fast <- as.numeric(fast)
  ref  <- as.numeric(ref)
  stopifnot(length(fast) == length(ref), length(fast) > 0)
  agree[[length(agree) + 1L]] <<- data.frame(
    method = method, reference = reference, quantity = quantity,
    n_compared = length(fast),
    max_abs_difference = max(abs(fast - ref)),
    max_abs_reference = max(abs(ref)),
    stringsAsFactors = FALSE
  )
}

# Kaplan-Meier survival at 1000 days.
ord  <- order(gbsg$rfstime)
fast <- survfit_fast(gbsg$rfstime[ord], gbsg$status[ord],
                     t_eval = 1000, conf.type = "log-log")
ref  <- summary(survfit(Surv(rfstime, status) ~ 1, data = gbsg),
                times = 1000)
add_agree("survfit_fast()", "survival::survfit()",
          "survival and standard error at 1000 days",
          unclass(fast)[c("surv", "std.err")], c(ref$surv, ref$std.err))

# Log-rank test, unstratified and stratified by tumor grade.
fast <- survdiff_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                      control = 0, side = 1)
ref  <- survdiff(Surv(rfstime, status) ~ hormon, data = gbsg)
add_agree("survdiff_fast()", "survival::survdiff()", "chi-square statistic",
          as.numeric(fast)^2, ref$chisq)

fast <- survdiff_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                      control = 0, side = 1, strata = gbsg$grade)
ref  <- survdiff(Surv(rfstime, status) ~ hormon + strata(grade), data = gbsg)
add_agree("survdiff_fast(strata)", "survival::survdiff()",
          "chi-square statistic, stratified by grade",
          as.numeric(fast)^2, ref$chisq)

# Fleming-Harrington weighted log-rank tests.
fh_grid <- data.frame(rho = c(0, 1, 0, 1), gamma = c(1, 0, 0, 1))
fast <- ref <- numeric(nrow(fh_grid))
for (i in seq_len(nrow(fh_grid))) {
  fast[i] <- as.numeric(survdiff_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                                      control = 0, side = 1, weight = "fh",
                                      rho = fh_grid$rho[i],
                                      gamma = fh_grid$gamma[i]))^2
  ref[i] <- nph::logrank.test(gbsg$rfstime, gbsg$status, gbsg$hormon,
                              rho = fh_grid$rho[i],
                              gamma = fh_grid$gamma[i])$test$Chisq
}
add_agree('survdiff_fast(weight = "fh")', "nph::logrank.test()",
          "chi-square statistic, FH(0,1), FH(1,0), FH(0,0), FH(1,1)",
          fast, ref)

# Cox hazard ratio (closed-form approximation of the Breslow partial
# likelihood), unstratified and stratified by tumor grade.
fast <- coxph_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                   control = 0, side = 1)
ref  <- coxph(Surv(rfstime, status) ~ hormon, data = gbsg, ties = "breslow")
add_agree("coxph_fast()", 'survival::coxph(ties = "breslow")',
          "log hazard ratio", unclass(fast)["coef"], coef(ref))

fast <- coxph_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                   control = 0, strata = gbsg$grade)
ref  <- coxph(Surv(rfstime, status) ~ hormon + strata(grade), data = gbsg,
              ties = "breslow")
add_agree("coxph_fast(strata)", 'survival::coxph(ties = "breslow")',
          "log hazard ratio, stratified by grade",
          unclass(fast)["coef"], coef(ref))

# Restricted mean survival time at 1000 days.
fast <- rmst_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                  control = 0, tau = 1000, side = 1)
ref  <- survRM2::rmst2(time = gbsg$rfstime, status = gbsg$status,
                       arm = gbsg$hormon, tau = 1000)
add_agree("rmst_fast()", "survRM2::rmst2()",
          "difference in RMST at 1000 days",
          unclass(fast)["diff"], ref$unadjusted.result[1, 1])

# Window mean survival time over [200, 1000] days against the integral of the
# Kaplan-Meier step function from survfit().
km_window_area <- function(gi, lo, hi) {
  sf <- survfit(Surv(rfstime, status) ~ 1, data = gbsg[gbsg$hormon == gi, ])
  tk <- c(0, sf$time)
  sk <- c(1, sf$surv)
  brk  <- sort(unique(c(lo, hi, sf$time[sf$time > lo & sf$time < hi])))
  area <- 0
  for (b in seq_len(length(brk) - 1L)) {
    u  <- brk[b]
    su <- sk[max(which(tk <= u))]
    area <- area + su * (brk[b + 1L] - u)
  }
  area
}
fast <- wmst_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                  control = 0, side = 1, tau1 = 200, tau2 = 1000)
add_agree("wmst_fast()", "Kaplan-Meier integral from survival::survfit()",
          "WMST per group over [200, 1000] days",
          unclass(fast)[c("wmst.control", "wmst.treatment")],
          c(km_window_area(0, 200, 1000), km_window_area(1, 200, 1000)))

# Weighted Kaplan-Meier test with the constant weight against the RMST
# difference over the observed range. rmst_fast() warns because tau extends
# past the shorter group's follow-up, which is intended here.
tmax <- max(gbsg$rfstime)
fast <- wkm_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                 control = 0, side = 1, weight = "constant")
ref  <- suppressWarnings(rmst_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                                   control = 0, tau = tmax, side = 1))
add_agree('wkm_fast(weight = "constant")', "rmst_fast() at the largest time",
          "weighted difference (identity with the RMST difference)",
          unclass(fast)["wdiff"], unclass(ref)["diff"])

# Milestone survival at 1000 days.
fast <- milestone_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                       control = 0, tau = 1000, method = "loglog", side = 1)
ref  <- summary(survfit(Surv(rfstime, status) ~ hormon, data = gbsg),
                times = 1000)
add_agree("milestone_fast()", "survival::survfit()",
          "survival per group at 1000 days",
          fast$surv[c("control", "treatment")], ref$surv)

# Kaplan-Meier median survival time.
fast <- medsurv_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                     control = 0, side = 1, method = "nph")
ref  <- summary(survfit(Surv(rfstime, status) ~ hormon,
                        data = gbsg))$table[, "median"]
add_agree("medsurv_fast()", "survival::survfit()", "median per group",
          unclass(fast)[c("median.control", "median.treatment")], ref)

# Average hazard with survival weight at 1000 days.
fast <- ahsw_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                  control = 0, tau = 1000, side = 1)
ref  <- survAH::ah2(time = gbsg$rfstime, status = gbsg$status,
                    arm = gbsg$hormon, tau = 1000)
add_agree("ahsw_fast()", "survAH::ah2()", "average hazard per group",
          unclass(fast)[c("ah.ctrl", "ah.trt")],
          c(ref$ah["AH (arm0)", "Est."], ref$ah["AH (arm1)", "Est."]))

# Average hazard ratio at 1000 days against a computation from the
# Kaplan-Meier curves of survfit().
fast <- ahr_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                 control = 0, tau = 1000, side = 1)
g    <- gbsg$hormon
ev   <- gbsg$rfstime[gbsg$status == 1]
grid <- sort(unique(c(0, ev[ev <= 1000], 1000)))
km_on_grid <- function(gi) {
  sf <- survfit(Surv(rfstime, status) ~ 1, data = gbsg[g == gi, ])
  approxfun(sf$time, sf$surv, method = "constant",
            yleft = 1, rule = 2, f = 0)(grid)
}
S0  <- km_on_grid(0)
S1  <- km_on_grid(1)
m   <- length(grid)
dS0 <- S0 - c(1, S0[-m])
GL  <- S0[m] * S1[m]
theta_ctrl <- -sum(S1 * dS0) / (1 - GL)
add_agree("ahr_fast()", "Kaplan-Meier curves from survival::survfit()",
          "group shares and average hazard ratio",
          c(fast$theta[[1]], fast$theta[[2]], fast$ahr),
          c(theta_ctrl, 1 - theta_ctrl, (1 - theta_ctrl) / theta_ctrl))

# Max-combo test: component Z statistics (the two packages orient the
# contrast in opposite directions, so absolute values are compared).
fast <- maxcombo_fast(gbsg$rfstime, gbsg$status, gbsg$hormon,
                      control = 0, side = 1, rho = c(0, 0, 1),
                      gamma = c(0, 1, 0))
ref  <- nph::logrank.maxtest(gbsg$rfstime, gbsg$status, gbsg$hormon)
add_agree("maxcombo_fast()", "nph::logrank.maxtest()",
          "absolute component Z statistics, FH(0,0), FH(0,1), FH(1,0)",
          abs(attr(fast, "z")), abs(ref$tests$z))

# Robust modestly-weighted log-rank test: component Z statistics.
gbsg_df <- data.frame(
  time  = gbsg$rfstime,
  event = gbsg$status,
  arm   = factor(ifelse(gbsg$hormon == 0, "control", "experimental"),
                 levels = c("control", "experimental"))
)
fast <- rmw_fast(gbsg_df$time, gbsg_df$event, gbsg_df$arm,
                 control = "control", side = 1, s_star = 0.5)
ref  <- c(nphRCT::wlrt(Surv(time, event) ~ arm, data = gbsg_df,
                       method = "mw", s_star = 1)$z,
          nphRCT::wlrt(Surv(time, event) ~ arm, data = gbsg_df,
                       method = "mw", s_star = 0.5)$z)
add_agree("rmw_fast()", "nphRCT::wlrt()",
          "absolute component Z statistics, s_star = 1 and 0.5",
          abs(attr(fast, "z")), abs(ref))

# Analysis cutoffs: "150 events, not before month 18 and not before 6 months
# after the 250th enrolled subject, and in any case by month 30", as in the
# validation vignette.
sim_c <- simdata_fast(nsim = 50, n = c(150, 150), a.time = c(0, 12),
                      a.rate = 300 / 12, e.median = list(12, 16),
                      d.hazard = 0.01, seed = 2026)
cut_c <- cutoff_fast(sim_c, event.looks = 150, time.looks = 18,
                     max.time = 30, min.enrolled = 250, min.followup = 6)
ref_c <- vapply(seq_len(50), function(s) {
  d <- sim_c[sim_c$sim == s, ]
  x <- data.frame(enroll_time = d$accrual_time,
                  cte = d$accrual_time + d$tte,
                  fail = d$event, stratum = "All")
  simtrial::get_analysis_date(x, planned_calendar_time = 18,
                              target_event_overall = 150,
                              max_extension_for_target_event = 30,
                              min_n_overall = 250, min_followup = 6)
}, numeric(1))
add_agree("cutoff_fast()", "simtrial::get_analysis_date()",
          "calendar time of the analysis in 50 simulated trials",
          cut_c[, 1], ref_c)

# analysis_fast() at two event-driven looks against survdiff() and rmst2()
# applied to the data censored at the same calendar cutoffs: a subject
# enrolled at a contributes if a <= cutoff, with observed time
# min(tte, cutoff - a) and an event if the event occurred by the cutoff.
res_a <- analysis_fast(sim_c, control = 1, event.looks = c(100, 180),
                       stat = c("logrank", "rmst"), tau = 6, side = 1)
res_a <- res_a[res_a$reached, ]
ref_a <- t(vapply(seq_len(nrow(res_a)), function(r) {
  d   <- sim_c[sim_c$sim == res_a$sim[r], ]
  cut <- res_a$cutoff[r]
  d   <- d[d$accrual_time <= cut, ]
  dd  <- data.frame(obs = pmin(d$tte, cut - d$accrual_time),
                    ev  = as.integer(d$event == 1 &
                                       d$accrual_time + d$tte <= cut),
                    grp = d$group)
  c(survdiff(Surv(obs, ev) ~ grp, data = dd)$chisq,
    survRM2::rmst2(time = dd$obs, status = dd$ev,
                   arm = as.integer(dd$grp == 2),
                   tau = 6)$unadjusted.result[1, 1])
}, numeric(2)))
add_agree("analysis_fast()", "survival::survdiff() on the censored data",
          "log-rank chi-square at two event-driven looks in 50 trials",
          res_a$logrank.chisq, ref_a[, 1])
add_agree("analysis_fast()", "survRM2::rmst2() on the censored data",
          "difference in RMST at two event-driven looks in 50 trials",
          res_a$rmst.diff, ref_a[, 2])

agreement <- do.call(rbind, agree)
rownames(agreement) <- NULL
print(agreement)
write.csv(agreement, file.path(out_dir, "bench_functions_agreement.csv"),
          row.names = FALSE)

# ---- Part 2: speed on a simulated data set ----------------------------------
dataset <- simdata_fast(nsim = 1, n = 500, a.time = c(0, 12.5), a.rate = 40,
                        e.median = list(5.811, 4.3), seed = 1)

# Sort once and reuse, the intended pattern for the pre-sorted fast path.
ord <- order(dataset$tte)
t_s <- dataset$tte[ord]
e_s <- dataset$event[ord]
g_s <- dataset$group[ord]

# Control is group 1, treatment is group 2.
arm <- as.integer(dataset$group == 2)

# Restriction horizon within both groups' follow-up, so survRM2 accepts it.
tau <- floor(min(tapply(t_s, g_s, max)))

df_rmw <- data.frame(
  tte   = dataset$tte,
  event = dataset$event,
  arm   = factor(ifelse(dataset$group == 1, "control", "treatment"),
                 levels = c("control", "treatment"))
)

B <- 1000

write.csv(data.frame(n = length(t_s), events = sum(e_s), tau = tau,
                     replicates = B),
          file.path(out_dir, "bench_functions_data.csv"), row.names = FALSE)

# Median times of the fast call and the reference, in milliseconds.
summarize_mb <- function(method, reference, mb) {
  s   <- summary(mb, unit = "ms")
  med <- stats::setNames(s$median, as.character(s$expr))
  data.frame(method = method, reference = reference,
             fast_ms = med[["fast"]], ref_ms = med[["ref"]],
             speedup = med[["ref"]] / med[["fast"]],
             stringsAsFactors = FALSE)
}
# Time two calls. The unevaluated call expressions are passed on with
# substitute(), so that microbenchmark evaluates them afresh in every replicate
# (passing them as ordinary arguments would time a cached promise).
mb <- function(fast, ref) {
  eval(substitute(microbenchmark::microbenchmark(fast = fast, ref = ref,
                                                 times = B)),
       parent.frame())
}

speed <- list(
  summarize_mb("survfit_fast()", "survival::survfit() + summary()", mb(
    survfit_fast(t_s, e_s, t_eval = tau, presorted = TRUE),
    summary(survfit(Surv(tte, event) ~ 1, data = dataset), times = tau))),
  summarize_mb("survdiff_fast()", "survival::survdiff()", mb(
    survdiff_fast(t_s, e_s, g_s, control = 1, side = 1, presorted = TRUE),
    survdiff(Surv(tte, event) ~ group, data = dataset))),
  summarize_mb("coxph_fast()", 'survival::coxph(ties = "breslow")', mb(
    coxph_fast(t_s, e_s, g_s, control = 1, side = 1, presorted = TRUE),
    coxph(Surv(tte, event) ~ I(group == 2), data = dataset,
          ties = "breslow"))),
  summarize_mb("rmst_fast()", "survRM2::rmst2()", mb(
    rmst_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
              presorted = TRUE),
    survRM2::rmst2(time = dataset$tte, status = dataset$event, arm = arm,
                   tau = tau))),
  summarize_mb('survdiff_fast(weight = "fh")', "nph::logrank.test()", mb(
    survdiff_fast(t_s, e_s, g_s, control = 1, side = 1, weight = "fh",
                  rho = 0, gamma = 1, presorted = TRUE),
    nph::logrank.test(dataset$tte, dataset$event, dataset$group,
                      rho = 0, gamma = 1))),
  summarize_mb("milestone_fast()", "survival::survfit() + summary()", mb(
    milestone_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
                   presorted = TRUE),
    summary(survfit(Surv(tte, event) ~ group, data = dataset), times = tau))),
  summarize_mb("medsurv_fast()", "nph::nphparams()", mb(
    medsurv_fast(t_s, e_s, g_s, control = 1, side = 1, method = "nph",
                 presorted = TRUE),
    nph::nphparams(time = dataset$tte, event = dataset$event, group = arm,
                   param_type = "Q", param_par = 0.5))),
  summarize_mb("maxcombo_fast()", "nph::logrank.maxtest()", mb(
    maxcombo_fast(t_s, e_s, g_s, control = 1, side = 1, rho = c(0, 0, 1),
                  gamma = c(0, 1, 0), presorted = TRUE),
    nph::logrank.maxtest(dataset$tte, dataset$event, arm))),
  summarize_mb("rmw_fast()", "nphRCT::wlrt() (two components)", mb(
    rmw_fast(t_s, e_s, g_s, control = 1, side = 1, s_star = 0.5,
             presorted = TRUE),
    {
      nphRCT::wlrt(Surv(tte, event) ~ arm, data = df_rmw, method = "mw",
                   s_star = 1)
      nphRCT::wlrt(Surv(tte, event) ~ arm, data = df_rmw, method = "mw",
                   s_star = 0.5)
    })),
  summarize_mb("ahsw_fast()", "survAH::ah2()", mb(
    ahsw_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
              presorted = TRUE),
    survAH::ah2(time = dataset$tte, status = dataset$event, arm = arm,
                tau = tau)))
)
speed <- do.call(rbind, speed)
rownames(speed) <- NULL
print(speed)
write.csv(speed, file.path(out_dir, "bench_functions_speed.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "bench_functions_sessionInfo.txt"))
