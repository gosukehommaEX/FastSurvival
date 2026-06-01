# Benchmark and profile the FastSurvival simulation trio
# (simdata_fast -> analysis_fast -> simsummary_fast).
#
# This is a local working script. It is not part of the package and does not
# need to pass R CMD check. Run it after devtools::load_all() in the package
# project root. It writes stage timings and an Rprof breakdown to ./profile_output/.
#
# Goal: determine where end-to-end time is spent so the fused-kernel design can
# target the real bottleneck (the hypothesis is that R-level orchestration in
# analysis_fast, plus data.frame materialisation in simdata_fast, dominates the
# fast C++ cores).

# ----------------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------------
if (!requireNamespace("FastSurvival", quietly = TRUE)) {
  stop("Load the package first, e.g. devtools::load_all().")
}
has_mb       <- requireNamespace("microbenchmark", quietly = TRUE)
has_simtrial <- requireNamespace("simtrial", quietly = TRUE)

out_dir <- file.path("tools", "profiling", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
log_path <- file.path(out_dir, "profile_log.txt")
cat("FastSurvival trio profile\n", file = log_path, append = FALSE)
cat(format(Sys.time()), "\n\n", file = log_path, append = TRUE)

log_line <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_path, append = TRUE)
}

# ----------------------------------------------------------------------------
# Scenario: two groups, NPH delayed effect, modest dropout, group-sequential
# looks. Control is exponential (median 12). Treatment shares the control
# hazard for the first 6 months, then drops (delayed separation). Enrolment is
# uniform over 12 time units; dropout is exponential (median 36).
# ----------------------------------------------------------------------------
build_design <- function(nsim, n_per_group = 200L) {
  list(
    nsim     = nsim,
    n        = c(n_per_group, n_per_group),
    a.time   = c(0, 12),
    a.rate   = 1,
    e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 18)),
    e.time   = list(NULL, c(0, 6, Inf)),
    d.median = list(36, 36),
    seed     = 2024L
  )
}

event_looks <- c(120, 180, 240)
stat_set    <- c("logrank", "coxph")
eff_bounds  <- c(-2.96, -2.36, -2.01)   # representative O'Brien-Fleming-like

run_generate <- function(nsim) {
  do.call(FastSurvival::simdata_fast, build_design(nsim))
}

run_analyze <- function(df) {
  FastSurvival::analysis_fast(df, control = 1, event.looks = event_looks,
                              stat = stat_set, side = 2)
}

run_summarize <- function(an) {
  FastSurvival::simsummary_fast(an, eff.col = "logrank.z",
                                efficacy = eff_bounds, direction = "lower")
}

# ----------------------------------------------------------------------------
# Stage timings across an nsim grid. Each stage is timed separately so the
# coarse split (generate vs analyze vs summarize) is visible. The median of a
# few repetitions is reported to be robust to noise.
# ----------------------------------------------------------------------------
time_median <- function(expr_fun, reps = 3L) {
  ts <- replicate(reps, {
    gc(FALSE)
    as.numeric(system.time(invisible(expr_fun()))[["elapsed"]])
  })
  median(ts)
}

nsim_grid <- c(500L, 1000L, 2000L, 5000L)

log_line("== Stage timings (median elapsed seconds) ==")
timing_rows <- vector("list", length(nsim_grid))
for (k in seq_along(nsim_grid)) {
  nsim <- nsim_grid[k]

  df <- NULL
  t_gen <- time_median(function() { df <<- run_generate(nsim) })

  an <- NULL
  t_ana <- time_median(function() { an <<- run_analyze(df) })

  t_sum <- time_median(function() { run_summarize(an) })

  t_tot <- t_gen + t_ana + t_sum
  n_rows <- nrow(df)

  timing_rows[[k]] <- data.frame(
    nsim      = nsim,
    data_rows = n_rows,
    generate  = round(t_gen, 4),
    analyze   = round(t_ana, 4),
    summarize = round(t_sum, 4),
    total     = round(t_tot, 4),
    analyze_pct = round(100 * t_ana / t_tot, 1)
  )
  log_line(sprintf(
    "nsim = %5d | rows = %8d | gen = %7.3f | ana = %7.3f | sum = %7.3f | tot = %7.3f | ana %% = %4.1f",
    nsim, n_rows, t_gen, t_ana, t_sum, t_tot, 100 * t_ana / t_tot))
}
timing_tab <- do.call(rbind, timing_rows)
write.csv(timing_tab, file.path(out_dir, "stage_timings.csv"), row.names = FALSE)
log_line("")

# ----------------------------------------------------------------------------
# Function-level profile at a representative nsim. Rprof attributes self time to
# split(), subsetting, data.frame construction, the per-look wrapper calls
# (survdiff_fast, coxph_fast, ...), and the C++ cores. This is the decisive
# diagnostic: if the cores are a small slice and the orchestration is large,
# the fused kernel is justified.
# ----------------------------------------------------------------------------
nsim_profile <- 4000L
log_line(sprintf("== Rprof breakdown at nsim = %d ==", nsim_profile))

prof_path <- file.path(out_dir, "trio.Rprof")
Rprof(prof_path, interval = 0.005, line.profiling = FALSE)
df_p <- run_generate(nsim_profile)
an_p <- run_analyze(df_p)
su_p <- run_summarize(an_p)
Rprof(NULL)

prof <- summaryRprof(prof_path)
self_tab  <- head(prof$by.self, 20L)
total_tab <- head(prof$by.total, 20L)

log_line(sprintf("Total profiled time: %.3f s (sampling interval 0.005 s)",
                 prof$sampling.time))
log_line("")
log_line("-- by.self (top 20) --")
self_txt <- capture.output(print(self_tab))
cat(self_txt, sep = "\n")
cat(self_txt, "\n", sep = "\n", file = log_path, append = TRUE)
log_line("")
log_line("-- by.total (top 20) --")
total_txt <- capture.output(print(total_tab))
cat(total_txt, sep = "\n")
cat(total_txt, "\n", sep = "\n", file = log_path, append = TRUE)
log_line("")

# Optional microbenchmark of the three stages at a single nsim for tighter
# variance estimates on the cheaper stages.
if (has_mb) {
  log_line("== microbenchmark (nsim = 2000, 10 iterations) ==")
  df_mb <- run_generate(2000L)
  an_mb <- run_analyze(df_mb)
  mb <- microbenchmark::microbenchmark(
    generate  = run_generate(2000L),
    analyze   = run_analyze(df_mb),
    summarize = run_summarize(an_mb),
    times     = 10L
  )
  mb_txt <- capture.output(print(mb))
  cat(mb_txt, sep = "\n")
  cat(mb_txt, "\n", sep = "\n", file = log_path, append = TRUE)
  log_line("")
}

# ----------------------------------------------------------------------------
# Optional simtrial end-to-end comparison. simtrial's API varies across
# versions, so this is wrapped defensively: any failure reports "skipped"
# rather than stopping the script. The comparison loops one trial per
# simulation (the simtrial idiom) to mirror a realistic OC loop.
# ----------------------------------------------------------------------------
bench_simtrial <- function(nsim) {
  enroll_rate  <- data.frame(duration = 12, rate = (2 * 200) / 12)
  fail_rate <- data.frame(
    stratum   = rep("All", 3),
    period    = c(1, 1, 2),
    treatment = c("control", "experimental", "experimental"),
    duration  = c(Inf, 6, Inf),
    rate      = c(log(2) / 12, log(2) / 12, log(2) / 18)
  )
  dropout_rate <- data.frame(
    stratum   = rep("All", 2),
    period    = c(1, 1),
    treatment = c("control", "experimental"),
    duration  = c(Inf, Inf),
    rate      = c(log(2) / 36, log(2) / 36)
  )

  t0 <- proc.time()[["elapsed"]]
  for (s in seq_len(nsim)) {
    dat <- simtrial::sim_pw_surv(
      n            = 2 * 200,
      enroll_rate  = enroll_rate,
      fail_rate    = fail_rate,
      dropout_rate = dropout_rate
    )
    cut <- simtrial::cut_data_by_event(dat, 240)
    invisible(simtrial::wlr(cut, weight = simtrial::fh(rho = 0, gamma = 0)))
  }
  proc.time()[["elapsed"]] - t0
}

if (has_simtrial) {
  log_line("== simtrial end-to-end (nsim = 2000), best effort ==")
  st_time <- tryCatch(bench_simtrial(2000L), error = function(e) {
    log_line(sprintf("simtrial comparison skipped: %s", conditionMessage(e)))
    NA_real_
  })
  if (is.finite(st_time)) {
    fs_time <- {
      df_c <- run_generate(2000L)
      g <- as.numeric(system.time(df_c <- run_generate(2000L))[["elapsed"]])
      a <- as.numeric(system.time(an_c <- run_analyze(df_c))[["elapsed"]])
      g + a
    }
    log_line(sprintf("simtrial (generate+cut+wlr loop): %.3f s", st_time))
    log_line(sprintf("FastSurvival (generate+analyze):  %.3f s", fs_time))
    log_line(sprintf("ratio simtrial / FastSurvival:    %.2fx", st_time / fs_time))
  }
} else {
  log_line("simtrial not installed: comparison skipped.")
}
log_line("")

# ----------------------------------------------------------------------------
# Record the execution environment for reproducibility.
# ----------------------------------------------------------------------------
log_line("== sessionInfo ==")
si_txt <- capture.output(print(sessionInfo()))
cat(si_txt, sep = "\n")
cat(si_txt, "\n", sep = "\n", file = log_path, append = TRUE)

log_line("")
log_line(sprintf("Outputs written to %s", normalizePath(out_dir)))
