# Crossover after a positive PFS analysis: FastSurvival (simdata_fast,
# cutoff_fast, analysis_fast, switch_fast) against TrialSimulator
# (CorrelatedPfsAndOs3 with a milestone crossover).
#
# Both simulate the same illness-death model with constant transition hazards.
# PFS is analyzed at d_pfs PFS events; when its one-sided log-rank p-value is at
# most alpha, control patients who are still on study switch at the later of
# their progression and the PFS analysis, and their remaining overall survival
# is multiplied by aft. OS is analyzed at d_os deaths. The script records the
# PFS power, the proportion of trials with crossover, the OS power, and the
# elapsed time, and writes them to tools/paper/output/bench_crossover.csv.
# The column distinct_trials counts the distinct calendar times of the PFS
# analysis, to confirm that the simulated trials are distinct.
#
# Run from the package root after installing FastSurvival 1.1.0 from CRAN
# (checked by machine_info.R):
#   source("tools/paper/data/bench_crossover.R")

library(FastSurvival)
source(file.path("tools", "paper", "data", "machine_info.R"))

out_dir <- file.path("tools", "paper", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Design ---------------------------------------------------------------
n_arm    <- 300
acc_dur  <- 18
acc_rate <- 2 * n_arm / acc_dur
h_ctrl   <- c(h01 = log(2) / 8,  h02 = log(2) / 30, h12 = log(2) / 12)
h_trt    <- c(h01 = log(2) / 12, h02 = log(2) / 36, h12 = log(2) / 12)
drop_haz <- -log(1 - 0.05) / 12
d_pfs    <- 380
d_os     <- 350
aft      <- 1.3
alpha    <- 0.025

nsim_fast <- 10000
nsim_ts   <- 2000

# One row of operating characteristics: the PFS and OS power and the mean
# calendar time of the two analyses, each with its Monte Carlo standard error
# over the simulated trials in which it is available, and the elapsed time.
prop_se <- function(x) {
  x <- x[!is.na(x)]
  sqrt(mean(x) * (1 - mean(x)) / length(x))
}
mean_se <- function(x) {
  x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}
oc_row <- function(package, nsim, pfs_rej, os_rej, pfs_cut, os_cut, elapsed) {
  data.frame(
    package = package, nsim = nsim,
    pfs_power = mean(pfs_rej, na.rm = TRUE), pfs_power_se = prop_se(pfs_rej),
    os_power = mean(os_rej, na.rm = TRUE), os_power_se = prop_se(os_rej),
    os_analyzed = sum(!is.na(os_rej)),
    pfs_month = mean(pfs_cut, na.rm = TRUE), pfs_month_se = mean_se(pfs_cut),
    os_month = mean(os_cut, na.rm = TRUE), os_month_se = mean_se(os_cut),
    elapsed = elapsed, sec_per_trial = elapsed / nsim,
    distinct_trials = length(unique(round(pfs_cut[!is.na(pfs_cut)], 8)))
  )
}

results <- list()

# ---- FastSurvival -----------------------------------------------------------
ep <- function(d, k) {
  data.frame(sim = d$sim, group = d$group, accrual_time = d$accrual_time,
             tte = d[[paste0("e", k, "_tte")]],
             event = d[[paste0("e", k, "_event")]])
}
t_fast <- system.time({
  df <- simdata_fast(nsim = nsim_fast, n = c(n_arm, n_arm),
                     a.time = c(0, acc_dur), a.rate = acc_rate,
                     h01.hazard = list(h_ctrl[["h01"]], h_trt[["h01"]]),
                     h02.hazard = list(h_ctrl[["h02"]], h_trt[["h02"]]),
                     h12.hazard = list(h_ctrl[["h12"]], h_trt[["h12"]]),
                     d.hazard = drop_haz, seed = 1)
  pc  <- cutoff_fast(df, event.looks = d_pfs, tte.col = "e1_tte",
                     event.col = "e1_event")
  pr  <- analysis_fast(ep(df, 1), control = 1, cutoff.looks = pc, side = 1)
  pos <- (pr$reached & pr$logrank.p <= alpha) %in% TRUE
  sw  <- switch_fast(df, group = 1, when = "later", cutoff = pc, sims = pos,
                     aft.factor = aft)
  oc  <- cutoff_fast(sw, event.looks = d_os, tte.col = "e2_tte",
                     event.col = "e2_event")
  orr <- analysis_fast(ep(sw, 2), control = 1, cutoff.looks = oc, side = 1)
})[["elapsed"]]
results$FastSurvival <- oc_row("FastSurvival", nsim_fast, pos,
                               orr$logrank.p <= alpha, pc[, 1], oc[, 1],
                               t_fast)

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
    mk_arm <- function(name, h) {
      a <- arm(name = name)
      a$add_endpoints(endpoint(name = c("pfs", "os"), type = c("tte", "tte"),
                               generator = CorrelatedPfsAndOs3,
                               h01 = h[["h01"]], h02 = h[["h02"]],
                               h12 = h[["h12"]]))
      a
    }
    tr <- trial(name = "crossover", n_patients = 2 * n_arm, duration = 500,
                enroller = StaggeredRecruiter,
                accrual_rate = data.frame(end_time = Inf,
                                          piecewise_rate = acc_rate),
                dropout = rexp, rate = drop_haz, seed = seed, silent = TRUE)
    tr$add_arms(sample_ratio = c(1, 1), mk_arm("control", h_ctrl),
                mk_arm("treatment", h_trt))

    what <- function(patient_data) {
      sw <- patient_data[patient_data$arm == "control", ]
      data.frame(patient_id = sw$patient_id, new_treatment = "treatment")
    }
    when <- function(patient_data) {
      data.frame(patient_id = patient_data$patient_id,
                 switch_time = pmax(patient_data$pfs,
                   patient_data$earliest_crossover_time_from_enrollment))
    }
    how <- function(patient_data) {
      data.frame(patient_id = patient_data$patient_id,
                 os = ifelse(patient_data$os > patient_data$switch_time,
                             patient_data$switch_time +
                               aft * (patient_data$os - patient_data$switch_time),
                             patient_data$os))
    }
    act_pfs <- function(trial) {
      d  <- trial$get_locked_data("pfs")
      lr <- fitLogrank(Surv(pfs, pfs_event) ~ arm, placebo = "control",
                       data = d, alternative = "less")
      trial$save(value = lr$p, name = "p_pfs")
      if (lr$p <= alpha) trial$crossover(what = what, how = how, when = when)
      invisible(NULL)
    }
    act_os <- function(trial) {
      d  <- trial$get_locked_data("os")
      lr <- fitLogrank(Surv(os, os_event) ~ arm, placebo = "control",
                       data = d, alternative = "less")
      trial$save(value = lr$p, name = "p_os")
      invisible(NULL)
    }
    lst <- listener(silent = TRUE)
    lst$add_milestones(
      milestone(name = "pfs", when = eventNumber(endpoint = "pfs", n = d_pfs),
                action = act_pfs),
      milestone(name = "os", when = eventNumber(endpoint = "os", n = d_os),
                action = act_os)
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
  results$TrialSimulator <- oc_row("TrialSimulator", nsim_ts,
                                   out$p_pfs <= alpha, out$p_os <= alpha,
                                   out[["milestone_time_<pfs>"]],
                                   out[["milestone_time_<os>"]], t_ts)
}

tab <- do.call(rbind, results)
rownames(tab) <- NULL
tab$speedup_vs_FastSurvival <- tab$sec_per_trial / tab$sec_per_trial[1]
# Differences from FastSurvival in units of their combined Monte Carlo standard
# error (0 in the FastSurvival row).
diff_z <- function(est, se) (est - est[1]) / sqrt(se^2 + se[1]^2)
tab$pfs_power_diff_z <- diff_z(tab$pfs_power, tab$pfs_power_se)
tab$os_power_diff_z  <- diff_z(tab$os_power, tab$os_power_se)
tab$pfs_month_diff_z <- diff_z(tab$pfs_month, tab$pfs_month_se)
tab$os_month_diff_z  <- diff_z(tab$os_month, tab$os_month_se)
print(tab)
write.csv(tab, file.path(out_dir, "bench_crossover.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "bench_crossover_sessionInfo.txt"))
