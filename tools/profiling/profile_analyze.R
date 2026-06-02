# Benchmark harness for the analyze stage of analysis_fast().
#
# Purpose: locate the dominant cost in the fused analysis kernel by measuring
# wall-clock time across simulation scales, look counts, and statistic sets.
# Because analysis_fast() is now a single .Call into analysis_loop_core, an
# R-level Rprof sees almost all self time inside .Call and cannot break the
# kernel down further. The look-count sweep is therefore the discriminator:
# the slope of time against the number of looks isolates the per-look
# administrative cut and time sort that v10 flagged as the remaining hot spot.
#
# Run from the package root in an interactive session:
#   source("tools/profiling/profile_analyze.R")
# Results are written under tools/profiling/output/.

suppressMessages(devtools::load_all("."))

out_dir <- file.path("tools", "profiling", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(out_dir, "analyze_timing_log.txt")
cat("", file = log_path)

log_msg <- function(...) {
  line <- paste0(..., "\n")
  cat(line)
  cat(line, file = log_path, append = TRUE)
}

# Median elapsed time over a few replications of a zero-argument function.
time_call <- function(f, reps = 5L) {
  ts <- numeric(reps)
  for (r in seq_len(reps)) {
    ts[r] <- system.time(f(), gcFirst = TRUE)[["elapsed"]]
  }
  stats::median(ts)
}

# Representative two-group data, optionally with a balanced subgroup factor.
make_data <- function(nsim, n_per, seed = 1L, with_subgroup = FALSE) {
  if (with_subgroup) {
    simdata_fast(nsim = nsim, n = c(n_per, n_per), a.time = c(0, 12),
                 a.rate = (2 * n_per) / 12,
                 e.hazard = list(list(0.10, 0.07), 0.05),
                 prevalence = c(0.5, 0.5), d.median = list(60, 60), seed = seed)
  } else {
    simdata_fast(nsim = nsim, n = c(n_per, n_per), a.time = c(0, 12),
                 a.rate = (2 * n_per) / 12,
                 e.hazard = list(0.07, 0.05),
                 d.median = list(60, 60), seed = seed)
  }
}

results <- data.frame()
record <- function(experiment, nsim, n_per, n_looks, n_pop, stat, look_kind,
                   seconds) {
  ncell <- nsim * n_looks * n_pop
  data.frame(experiment = experiment, nsim = nsim, n_per = n_per,
             n_looks = n_looks, n_pop = n_pop, stat = stat,
             look_kind = look_kind, seconds = round(seconds, 4),
             us_per_cell = round(1e6 * seconds / ncell, 3),
             stringsAsFactors = FALSE)
}

log_msg("== analysis_fast analyze-stage benchmark ==")
log_msg("Started: ", format(Sys.time()))

# --- Experiment 1: look-count sweep (per-look cost discriminator) ----------
log_msg("\n[1] Look-count sweep (nsim = 5000, n_per = 175, stat = logrank)")
dat1 <- make_data(nsim = 5000, n_per = 175, seed = 11)
for (L in c(1L, 2L, 3L, 5L)) {
  looks_t <- seq(12, 36, length.out = L)
  tt <- time_call(function() {
    analysis_fast(dat1, control = 1, time.looks = looks_t, stat = "logrank")
  })
  results <- rbind(results,
    record("looks_sweep", 5000, 175, L, 1L, "logrank", "time", tt))
  log_msg(sprintf("  looks = %d : %.4f s (%.3f us/cell)",
                  L, tt, 1e6 * tt / (5000 * L)))
}

# --- Experiment 2: statistic sets at a fixed two-look configuration --------
log_msg("\n[2] Statistic sets (nsim = 5000, n_per = 175, time.looks = 18, 30)")
stat_sets <- list(
  logrank     = "logrank",
  coxph       = "coxph",
  rmst        = "rmst",
  km          = "km",
  lr_cox_rmst = c("logrank", "coxph", "rmst"),
  maxcombo    = "maxcombo",
  ahsw        = "ahsw"
)
for (nm in names(stat_sets)) {
  s     <- stat_sets[[nm]]
  tau_v <- if (any(c("rmst", "ahsw") %in% s)) 18 else NULL
  te_v  <- if ("km" %in% s) 18 else NULL
  reps  <- if ("maxcombo" %in% s) 3L else 5L
  tt <- time_call(function() {
    analysis_fast(dat1, control = 1, time.looks = c(18, 30), stat = s,
                  tau = tau_v, t.eval = te_v)
  }, reps = reps)
  results <- rbind(results,
    record("stat_sets", 5000, 175, 2L, 1L, nm, "time", tt))
  log_msg(sprintf("  %-12s : %.4f s (%.3f us/cell)",
                  nm, tt, 1e6 * tt / (5000 * 2)))
}

# --- Experiment 3: scale sweep over nsim and n_per -------------------------
log_msg("\n[3] Scale sweep (stat = logrank, time.looks = 18, 30)")
for (np in c(100L, 350L)) {
  for (ns in c(2000L, 5000L, 10000L)) {
    d <- make_data(nsim = ns, n_per = np, seed = 21)
    tt <- time_call(function() {
      analysis_fast(d, control = 1, time.looks = c(18, 30), stat = "logrank")
    })
    results <- rbind(results,
      record("scale_sweep", ns, np, 2L, 1L, "logrank", "time", tt))
    log_msg(sprintf("  nsim = %5d, n_per = %3d : %.4f s (%.3f us/cell)",
                    ns, np, tt, 1e6 * tt / (ns * 2)))
  }
}

# --- Experiment 4: event-driven versus calendar looks ----------------------
log_msg("\n[4] Event-driven vs calendar looks (nsim = 5000, n_per = 175)")
tt_time <- time_call(function() {
  analysis_fast(dat1, control = 1, time.looks = c(18, 24, 30),
                stat = "logrank")
})
tt_event <- time_call(function() {
  analysis_fast(dat1, control = 1, event.looks = c(120, 180, 240),
                stat = "logrank")
})
results <- rbind(results,
  record("look_kind", 5000, 175, 3L, 1L, "logrank", "time", tt_time),
  record("look_kind", 5000, 175, 3L, 1L, "logrank", "event", tt_event))
log_msg(sprintf("  calendar : %.4f s", tt_time))
log_msg(sprintf("  event    : %.4f s", tt_event))

# --- Experiment 5: subgroup marginalization on or off ----------------------
log_msg("\n[5] Subgroup off vs on (nsim = 5000, n_per = 175, two looks)")
dat5 <- make_data(nsim = 5000, n_per = 175, seed = 31, with_subgroup = TRUE)
tt_off <- time_call(function() {
  analysis_fast(dat5, control = 1, time.looks = c(18, 30), stat = "logrank")
})
tt_on <- time_call(function() {
  analysis_fast(dat5, control = 1, time.looks = c(18, 30), stat = "logrank",
                by.subgroup = TRUE)
})
n_pop_on <- 1L + length(unique(dat5$subgroup))
results <- rbind(results,
  record("subgroup", 5000, 175, 2L, 1L,        "logrank", "time", tt_off),
  record("subgroup", 5000, 175, 2L, n_pop_on,  "logrank", "time", tt_on))
log_msg(sprintf("  by.subgroup = FALSE : %.4f s", tt_off))
log_msg(sprintf("  by.subgroup = TRUE  : %.4f s (n_pop = %d)", tt_on, n_pop_on))

# --- Experiment 6: R-level Rprof sanity check (expected uninformative) ------
log_msg("\n[6] R-level Rprof of one analyze call (self time concentrates in .Call)")
prof_path <- file.path(out_dir, "analyze_rprof.out")
Rprof(prof_path, interval = 0.005)
invisible(analysis_fast(dat1, control = 1, time.looks = c(18, 24, 30),
                        stat = "logrank"))
Rprof(NULL)
prof_summary <- summaryRprof(prof_path)
log_msg("  Top self-time entries (by.self):")
top_self <- utils::head(prof_summary$by.self, 6)
for (i in seq_len(nrow(top_self))) {
  log_msg(sprintf("    %-24s self = %5.1f%%",
                  rownames(top_self)[i], top_self$self.pct[i]))
}

# --- Persist results -------------------------------------------------------
saveRDS(results, file = file.path(out_dir, "analyze_timing.rds"))
log_msg("\n== Results table ==")
log_msg(paste(utils::capture.output(print(results, row.names = FALSE)),
              collapse = "\n"))

log_msg("\n== sessionInfo() ==")
log_msg(paste(utils::capture.output(sessionInfo()), collapse = "\n"))
log_msg("\nFinished: ", format(Sys.time()))
