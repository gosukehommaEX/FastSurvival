# End-to-end comparison of FastSurvival, simtrial, and TrialSimulator on the
# same two-arm group-sequential design (log-rank test, two event-driven looks).
#
# For each package the script records the operating characteristics (power,
# expected number of events and calendar time at each look) and the elapsed
# time, and writes them to tools/paper/output/bench_gsd.csv. The operating
# characteristics should agree within Monte Carlo error; the elapsed time per
# simulated trial is the benchmark. The column distinct_trials counts the
# distinct calendar times of the final look, to confirm that the simulated
# trials are distinct.
#
# Run from the package root (FastSurvival.Rproj) after installing FastSurvival
# 1.1.0 from CRAN (checked by machine_info.R):
#   source("tools/paper/data/bench_gsd.R")

library(FastSurvival)
source(file.path("tools", "paper", "data", "machine_info.R"))

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

# One row of operating characteristics: the power and the mean Z statistic and
# calendar time of each look, each with its Monte Carlo standard error, and the
# elapsed time.
mean_se <- function(x) sd(x) / sqrt(length(x))
oc_row <- function(package, z, cut1, cut2, elapsed) {
  rej <- reject_gs(z, bnd)
  power <- mean(rej)
  data.frame(
    package = package, nsim = nrow(z),
    power = power, power_se = sqrt(power * (1 - power) / length(rej)),
    mean_z_look1 = mean(z[, 1]), mean_z_look1_se = mean_se(z[, 1]),
    mean_z_look2 = mean(z[, 2]), mean_z_look2_se = mean_se(z[, 2]),
    cutoff_look1 = mean(cut1), cutoff_look1_se = mean_se(cut1),
    cutoff_look2 = mean(cut2), cutoff_look2_se = mean_se(cut2),
    elapsed = elapsed, sec_per_trial = elapsed / nrow(z),
    distinct_trials = length(unique(round(cut2, 8)))
  )
}

nsim_fast <- 10000
nsim_simtrial <- 2000
nsim_ts <- 2000

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
results$FastSurvival <- oc_row("FastSurvival", z_fast,
                               res$cutoff[res$look == 1],
                               res$cutoff[res$look == 2], t_fast)

# ---- simtrial ---------------------------------------------------------------
if (requireNamespace("simtrial", quietly = TRUE)) {
  # Run simtrial sequentially, as the other packages, and reproducibly.
  if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
  set.seed(1)
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
  results$simtrial <- oc_row("simtrial", z_st,
                             st$cut_date[st$analysis == 1],
                             st$cut_date[st$analysis == 2], t_st)
}

# ---- TrialSimulator ---------------------------------------------------------
if (requireNamespace("TrialSimulator", quietly = TRUE)) {
  suppressPackageStartupMessages(library(TrialSimulator))
  # TrialSimulator 1.35.8 draws the seed of each replicate of
  # controller$run(n) from the random-number state left by the previous
  # replicate, so the seeds follow a deterministic sequence that can repeat;
  # with seed = 1 the replicates repeated after a few hundred trials
  # (check_crossover_pfs.R). Each simulated trial is therefore run by its own
  # controller with seed = i, and only the time of controller$run() is
  # counted.
  make_ctl <- function(seed) {
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
                dropout = rexp, rate = drop_haz, seed = seed, silent = TRUE)
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
    controller(tr, lst)
  }
  t_ts <- 0
  outs <- vector("list", nsim_ts)
  for (i in seq_len(nsim_ts)) {
    ctl <- make_ctl(i)
    t_ts <- t_ts + system.time(
      ctl$run(n = 1, plot_event = FALSE, silent = TRUE))[["elapsed"]]
    outs[[i]] <- ctl$get_output()
  }
  out <- as.data.frame(dplyr::bind_rows(outs))
  # fitLogrank's z has the sign of the log hazard ratio (negative for benefit).
  z_ts <- -cbind(out$z_ia, out$z_fa)
  results$TrialSimulator <- oc_row("TrialSimulator", z_ts,
                                   out[["milestone_time_<ia>"]],
                                   out[["milestone_time_<fa>"]], t_ts)
}

tab <- do.call(rbind, results)
rownames(tab) <- NULL
tab$speedup_vs_FastSurvival <- tab$sec_per_trial / tab$sec_per_trial[1]
# Differences from FastSurvival in units of their combined Monte Carlo standard
# error (0 in the FastSurvival row).
diff_z <- function(est, se) (est - est[1]) / sqrt(se^2 + se[1]^2)
tab$power_diff_z   <- diff_z(tab$power, tab$power_se)
tab$z_look1_diff_z <- diff_z(tab$mean_z_look1, tab$mean_z_look1_se)
tab$z_look2_diff_z <- diff_z(tab$mean_z_look2, tab$mean_z_look2_se)
tab$cutoff1_diff_z <- diff_z(tab$cutoff_look1, tab$cutoff_look1_se)
tab$cutoff2_diff_z <- diff_z(tab$cutoff_look2, tab$cutoff_look2_se)
print(tab)
write.csv(tab, file.path(out_dir, "bench_gsd.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "bench_gsd_sessionInfo.txt"))
