# End-to-end comparison of FastSurvival, simtrial, and TrialSimulator on the
# same two-arm group-sequential design (log-rank test, two event-driven looks).
#
# For each package the script records the operating characteristics (power,
# expected number of events and calendar time at each look) and the elapsed
# time, and writes them to tools/paper/output/bench_gsd.csv. The operating
# characteristics should agree within Monte Carlo error; the elapsed time per
# simulated trial is the benchmark.
#
# Run from the package root (FastSurvival.Rproj) after installing the package:
#   source("tools/paper/bench_gsd.R")

library(FastSurvival)

out_dir <- file.path("tools", "paper", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Design ---------------------------------------------------------------
n_arm      <- 300                       # per arm
acc_dur    <- 12                        # months
acc_rate   <- 2 * n_arm / acc_dur       # subjects per month
med_ctrl   <- 12
hr         <- 0.75
drop_haz   <- 0.001                     # per month
looks      <- c(200, 350)               # cumulative events
alpha      <- 0.025

bnd <- gsDesign::gsDesign(k = 2, test.type = 1, alpha = alpha,
                          timing = looks / max(looks),
                          sfu = gsDesign::sfLDOF)$upper$bound

# Sequential rejection with efficacy bounds on the benefit-positive Z scale.
reject_gs <- function(zmat, bound) {
  rej <- rep(FALSE, nrow(zmat))
  for (k in seq_len(ncol(zmat))) rej <- rej | (zmat[, k] >= bound[k])
  rej
}

nsim_fast <- 10000
nsim_simtrial <- 1000
nsim_ts <- 200

results <- list()

# ---- FastSurvival -----------------------------------------------------------
t_fast <- system.time({
  df <- simdata_fast(nsim = nsim_fast, n = c(n_arm, n_arm),
                     a.time = c(0, acc_dur), a.rate = acc_rate,
                     e.median = list(med_ctrl, med_ctrl / hr),
                     d.hazard = drop_haz, seed = 1)
  res <- analysis_fast(df, control = 1, event.looks = looks, side = 1)
})[["elapsed"]]
# logrank.z is negative for benefit; flip to the benefit-positive scale.
z_fast <- -matrix(res$logrank.z, ncol = length(looks), byrow = TRUE)
results$FastSurvival <- data.frame(
  package = "FastSurvival", nsim = nsim_fast,
  power = mean(reject_gs(z_fast, bnd)),
  mean_z_look1 = mean(z_fast[, 1]), mean_z_look2 = mean(z_fast[, 2]),
  cutoff_look1 = mean(res$cutoff[res$look == 1]),
  cutoff_look2 = mean(res$cutoff[res$look == 2]),
  elapsed = t_fast, sec_per_trial = t_fast / nsim_fast
)

# ---- simtrial ---------------------------------------------------------------
if (requireNamespace("simtrial", quietly = TRUE)) {
  t_st <- system.time({
    st <- simtrial::sim_gs_n(
      n_sim = nsim_simtrial, sample_size = 2 * n_arm,
      enroll_rate = data.frame(duration = acc_dur, rate = acc_rate),
      fail_rate = data.frame(stratum = "All", duration = 1000,
                             fail_rate = log(2) / med_ctrl, hr = hr,
                             dropout_rate = drop_haz),
      test = simtrial::wlr,
      cut = list(ia = simtrial::create_cut(target_event_overall = looks[1]),
                 fa = simtrial::create_cut(target_event_overall = looks[2])),
      weight = simtrial::fh(rho = 0, gamma = 0)
    )
  })[["elapsed"]]
  st <- st[order(st$sim_id, st$analysis), ]
  z_st <- matrix(st$z, ncol = length(looks), byrow = TRUE)
  # simtrial reports z = -estimate / se; check that benefit is positive.
  if (mean(z_st[, 2]) < 0) z_st <- -z_st
  results$simtrial <- data.frame(
    package = "simtrial", nsim = nsim_simtrial,
    power = mean(reject_gs(z_st, bnd)),
    mean_z_look1 = mean(z_st[, 1]), mean_z_look2 = mean(z_st[, 2]),
    cutoff_look1 = mean(st$cut_date[st$analysis == 1]),
    cutoff_look2 = mean(st$cut_date[st$analysis == 2]),
    elapsed = t_st, sec_per_trial = t_st / nsim_simtrial
  )
}

# ---- TrialSimulator ---------------------------------------------------------
if (requireNamespace("TrialSimulator", quietly = TRUE)) {
  suppressPackageStartupMessages(library(TrialSimulator))
  t_ts <- system.time({
    ctrl <- arm(name = "control")
    ctrl$add_endpoints(endpoint(name = "os", type = "tte", generator = rexp,
                                rate = log(2) / med_ctrl))
    trt <- arm(name = "treatment")
    trt$add_endpoints(endpoint(name = "os", type = "tte", generator = rexp,
                               rate = log(2) / med_ctrl * hr))
    tr <- trial(name = "gsd", n_patients = 2 * n_arm, duration = 500,
                enroller = StaggeredRecruiter,
                accrual_rate = data.frame(end_time = Inf,
                                          piecewise_rate = acc_rate),
                dropout = rexp, rate = drop_haz, seed = 1, silent = TRUE)
    tr$add_arms(sample_ratio = c(1, 1), ctrl, trt)
    act <- function(look) {
      force(look)
      function(trial) {
        d  <- trial$get_locked_data(look)
        lr <- fitLogrank(Surv(os, os_event) ~ arm, placebo = "control",
                         data = d, alternative = "less")
        trial$save(value = lr$z, name = paste0("z_", look))
        invisible(NULL)
      }
    }
    lst <- listener(silent = TRUE)
    lst$add_milestones(
      milestone(name = "ia", when = eventNumber(endpoint = "os", n = looks[1]),
                action = act("ia")),
      milestone(name = "fa", when = eventNumber(endpoint = "os", n = looks[2]),
                action = act("fa"))
    )
    ctl <- controller(tr, lst)
    ctl$run(n = nsim_ts, plot_event = FALSE, silent = TRUE)
    out <- ctl$get_output()
  })[["elapsed"]]
  # fitLogrank's z has the sign of the log hazard ratio (negative for benefit).
  z_ts <- -cbind(out$z_ia, out$z_fa)
  results$TrialSimulator <- data.frame(
    package = "TrialSimulator", nsim = nsim_ts,
    power = mean(reject_gs(z_ts, bnd)),
    mean_z_look1 = mean(z_ts[, 1]), mean_z_look2 = mean(z_ts[, 2]),
    cutoff_look1 = mean(out[["milestone_time_<ia>"]]),
    cutoff_look2 = mean(out[["milestone_time_<fa>"]]),
    elapsed = t_ts, sec_per_trial = t_ts / nsim_ts
  )
}

tab <- do.call(rbind, results)
rownames(tab) <- NULL
tab$speedup_vs_FastSurvival <- tab$sec_per_trial / tab$sec_per_trial[1]
print(tab)
write.csv(tab, file.path(out_dir, "bench_gsd.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "bench_gsd_sessionInfo.txt"))
