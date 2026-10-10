# Check of the PFS analysis of bench_crossover.R.
#
# In bench_crossover.R the PFS power of TrialSimulator (2,000 trials) exceeded
# that of FastSurvival (10,000 trials) by about four Monte Carlo standard
# errors, while an independent simulation in Python agreed with FastSurvival.
# This script records, for every trial simulated by TrialSimulator with and
# without the crossover action, the calendar time of the PFS analysis, the
# number of PFS events and of subjects per group in the analyzed data, the
# log-rank Z statistic, and the crude log hazard ratio (log of the ratio of
# events per unit of follow-up, the maximum likelihood estimate under
# exponential PFS). It also records the same quantities for FastSurvival from
# 100,000 trials generated in 10 batches from independent dqrng streams.
# Results are written to tools/paper/output/check_crossover_pfs_ts.csv and
# check_crossover_pfs_fs.csv.
#
# Run from the package root after installing FastSurvival 1.2.0 from CRAN
# (checked by machine_info.R):
#   source("tools/paper/data/check_crossover_pfs.R")
# or from the article folder with source("scripts/check_crossover_pfs.R").

library(FastSurvival)
# The scripts are in tools/paper/data of the package or in scripts of the
# article folder; machine_info.R sets the output folder paper_out_dir.
paper_script_dir <- if (dir.exists(file.path("tools", "paper", "data"))) {
  file.path("tools", "paper", "data")
} else if (file.exists(file.path("scripts", "machine_info.R"))) {
  "scripts"
} else {
  stop("Run the script from the package root or from the article folder.",
       call. = FALSE)
}
source(file.path(paper_script_dir, "machine_info.R"))

out_dir <- paper_out_dir

# ---- Design (as in bench_crossover.R) ---------------------------------------
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

# ---- FastSurvival: 10 batches of 10,000 trials ------------------------------
ep <- function(d, k) {
  data.frame(sim = d$sim, group = d$group, accrual_time = d$accrual_time,
             tte = d[[paste0("e", k, "_tte")]],
             event = d[[paste0("e", k, "_event")]])
}
fs_batches <- paper_n(10, 2)
fs <- do.call(rbind, lapply(seq_len(fs_batches), function(b) {
  df <- simdata_fast(nsim = paper_n(10000, 200), n = c(n_arm, n_arm),
                     a.time = c(0, acc_dur), a.rate = acc_rate,
                     h01.hazard = list(h_ctrl[["h01"]], h_trt[["h01"]]),
                     h02.hazard = list(h_ctrl[["h02"]], h_trt[["h02"]]),
                     h12.hazard = list(h_ctrl[["h12"]], h_trt[["h12"]]),
                     d.hazard = drop_haz, seed = 2, stream = b)
  pc <- cutoff_fast(df, event.looks = d_pfs, tte.col = "e1_tte",
                    event.col = "e1_event")
  pf <- ep(df, 1)
  pr <- analysis_fast(pf, control = 1, cutoff.looks = pc, side = 1)
  # Crude log hazard ratio of the data censored at the PFS cutoff.
  cut <- pc[match(pf$sim, sort(unique(pf$sim))), 1]
  inc <- !is.na(cut) & pf$accrual_time <= cut
  obs <- pmin(pf$tte, cut - pf$accrual_time)[inc]
  evt <- (pf$event == 1 & pf$accrual_time + pf$tte <= cut)[inc]
  key <- paste(pf$sim[inc], pf$group[inc])
  ev  <- rowsum(as.numeric(evt), key)
  pt  <- rowsum(obs, key)
  sim_g <- do.call(rbind, strsplit(rownames(ev), " "))
  rate  <- ev[, 1] / pt[, 1]
  lhr <- log(rate[sim_g[, 2] == "2"] / rate[sim_g[, 2] == "1"])
  rej <- (pr$reached & pr$logrank.p <= alpha) %in% TRUE
  rm(df, pf)
  gc()
  data.frame(batch = b, nsim = nrow(pr), rejections = sum(rej),
             sum_z = sum(pr$logrank.z), sum_z2 = sum(pr$logrank.z^2),
             sum_cutoff = sum(pc[, 1]), sum_cutoff2 = sum(pc[, 1]^2),
             min_events = min(pr$n.event), max_events = max(pr$n.event),
             sum_log_hr = sum(lhr), sum_log_hr2 = sum(lhr^2))
}))
print(fs)
write.csv(fs, file.path(out_dir, "check_crossover_pfs_fs.csv"),
          row.names = FALSE)

# ---- TrialSimulator: 2,000 trials with and without the crossover ------------
if (requireNamespace("TrialSimulator", quietly = TRUE)) {
  suppressPackageStartupMessages(library(TrialSimulator))
  run_ts <- function(crossover, nsim) {
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
                dropout = rexp, rate = drop_haz, seed = 1, silent = TRUE)
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
                               aft * (patient_data$os -
                                        patient_data$switch_time),
                             patient_data$os))
    }
    act_pfs <- function(trial) {
      d  <- trial$get_locked_data("pfs")
      lr <- fitLogrank(Surv(pfs, pfs_event) ~ arm, placebo = "control",
                       data = d, alternative = "less", tidy = FALSE)
      ev <- tapply(d$pfs_event, d$arm, sum)
      pt <- tapply(d$pfs, d$arm, sum)
      trial$save(value = lr$p, name = "chk_p")
      trial$save(value = lr$z, name = "chk_z")
      trial$save(value = lr$info, name = "chk_events")
      trial$save(value = lr$n_pbo, name = "chk_n_ctrl")
      trial$save(value = lr$n_trt, name = "chk_n_trt")
      trial$save(value = nrow(d), name = "chk_rows")
      trial$save(value = log((ev[["treatment"]] / pt[["treatment"]]) /
                               (ev[["control"]] / pt[["control"]])),
                 name = "chk_log_hr")
      if (crossover && lr$p <= alpha) {
        trial$crossover(what = what, how = how, when = when)
      }
      invisible(NULL)
    }
    act_os <- function(trial) {
      d  <- trial$get_locked_data("os")
      lr <- fitLogrank(Surv(os, os_event) ~ arm, placebo = "control",
                       data = d, alternative = "less")
      trial$save(value = lr$p, name = "chk_p_os")
      invisible(NULL)
    }
    lst <- listener(silent = TRUE)
    lst$add_milestones(
      milestone(name = "pfs", when = eventNumber(endpoint = "pfs", n = d_pfs),
                action = act_pfs),
      milestone(name = "os", when = eventNumber(endpoint = "os", n = d_os),
                action = act_os)
    )
    ctl <- controller(tr, lst)
    ctl$run(n = nsim, plot_event = FALSE, silent = TRUE)
    out <- ctl$get_output()
    data.frame(crossover = crossover, trial = seq_len(nrow(out)),
               cutoff = out[["milestone_time_<pfs>"]],
               events = out$chk_events, rows = out$chk_rows,
               n_ctrl = out$chk_n_ctrl, n_trt = out$chk_n_trt,
               z = out$chk_z, p = out$chk_p, log_hr = out$chk_log_hr,
               p_os = out$chk_p_os)
  }
  n_ts <- paper_n(2000, 50)
  ts <- rbind(run_ts(FALSE, n_ts), run_ts(TRUE, n_ts))
  write.csv(ts, file.path(out_dir, "check_crossover_pfs_ts.csv"),
            row.names = FALSE)
  print(tapply(ts$p <= alpha, ts$crossover, mean))
  print(aggregate(ts[, c("z", "events", "n_ctrl", "n_trt", "log_hr",
                         "cutoff")],
                  list(crossover = ts$crossover), mean))
}
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "check_crossover_pfs_sessionInfo.txt"))
