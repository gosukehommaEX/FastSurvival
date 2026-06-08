# Speed benchmark for FastSurvival against established reference packages.
#
# Generates a single two-group dataset with simdata_fast(), then times each of
# the twelve estimation and testing functions against an established reference
# with microbenchmark (1000 replicates each). FastSurvival calls use one-sided
# tests (side = 1) where applicable and the pre-sorted fast path
# (presorted = TRUE), reflecting use inside a simulation loop where the data are
# sorted once and reused. References on CRAN are survival, survRM2, survAH, nph
# and nphRCT; survWMST, nphsim and AHR are on GitHub, and those rows are skipped
# when the package is not installed. Results are printed and written to
# speed_results.csv / speed_results.rds under tools/benchmarks/.
#
# Run from the package root, for example: Rscript tools/benchmark_speed.R

library(FastSurvival)

out_dir <- file.path("tools", "benchmarks")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

log_path <- file.path(out_dir, "speed_results.log")
log_line <- function(...) {
  cat(paste0(..., "\n"), file = log_path, append = TRUE)
}
cat("FastSurvival speed benchmark\n", file = log_path, append = FALSE)
log_line("Run started: ", format(Sys.time()))

have_survival <- requireNamespace("survival", quietly = TRUE)
have_survRM2  <- requireNamespace("survRM2", quietly = TRUE)
have_survAH   <- requireNamespace("survAH", quietly = TRUE)
have_nph      <- requireNamespace("nph", quietly = TRUE)
have_nphRCT   <- requireNamespace("nphRCT", quietly = TRUE)
have_survWMST <- requireNamespace("survWMST", quietly = TRUE)
have_nphsim   <- requireNamespace("nphsim", quietly = TRUE)
have_AHR      <- requireNamespace("AHR", quietly = TRUE)
have_mb       <- requireNamespace("microbenchmark", quietly = TRUE)

if (!have_mb) {
  stop("Package 'microbenchmark' is required to run this benchmark.")
}
if (have_survival) {
  library(survival)
}

# ---- Benchmark dataset -------------------------------------------------------
dataset <- simdata_fast(
  nsim     = 1,
  n        = 500,
  a.time   = c(0, 12.5),
  a.rate   = 40,
  e.median = list(5.811, 4.3),
  seed     = 1
)

# Sort once and reuse, the intended pattern for the pre-sorted fast path.
ord <- order(dataset$tte)
t_s <- dataset$tte[ord]
e_s <- dataset$event[ord]
g_s <- dataset$group[ord]

# Control is group 1, treatment is group 2.
arm <- as.integer(dataset$group == 2)

# Group coded as 0/1 for nph functions (which require 0/1).
group0 <- dataset$group - 1

# Restriction horizon within both arms' follow-up, so survRM2 accepts it.
tau <- floor(min(tapply(t_s, g_s, max)))

# Factor arm for the nphRCT reference used in the rmw_fast benchmark.
df_rmw <- data.frame(
  tte   = dataset$tte,
  event = dataset$event,
  arm   = factor(ifelse(dataset$group == 1, "control", "treatment"),
                 levels = c("control", "treatment"))
)

B <- 1000

log_line("n = ", length(t_s), ", events = ", sum(e_s), ", tau = ", tau,
         ", replicates = ", B)

# ---- Helper: median fast and reference times, in milliseconds ----------------
summarize_mb <- function(label, replaces, mb) {
  s   <- summary(mb, unit = "ms")
  med <- stats::setNames(s$median, as.character(s$expr))
  data.frame(
    method   = label,
    replaces = replaces,
    fast_ms  = round(med[["fast"]], 4),
    ref_ms   = round(med[["ref"]], 4),
    speedup  = round(med[["ref"]] / med[["fast"]], 1),
    stringsAsFactors = FALSE
  )
}

results <- list()

# ---- survfit_fast vs survfit + summary --------------------------------------
if (have_survival) {
  results[["survfit"]] <- summarize_mb(
    "survfit_fast()", "survfit() + summary()",
    microbenchmark::microbenchmark(
      fast = survfit_fast(t_s, e_s, t_eval = tau, presorted = TRUE),
      ref  = summary(survfit(Surv(tte, event) ~ 1, data = dataset),
                     times = tau),
      times = B
    )
  )
}

# ---- survdiff_fast vs survdiff ----------------------------------------------
if (have_survival) {
  results[["survdiff"]] <- summarize_mb(
    "survdiff_fast()", "survdiff()",
    microbenchmark::microbenchmark(
      fast = survdiff_fast(t_s, e_s, g_s, control = 1, side = 1,
                           presorted = TRUE),
      ref  = survdiff(Surv(tte, event) ~ group, data = dataset),
      times = B
    )
  )
}

# ---- coxph_fast vs coxph -----------------------------------------------------
if (have_survival) {
  results[["coxph"]] <- summarize_mb(
    "coxph_fast()", "coxph()",
    microbenchmark::microbenchmark(
      fast = coxph_fast(t_s, e_s, g_s, control = 1, side = 1,
                        presorted = TRUE),
      ref  = coxph(Surv(tte, event) ~ I(group == 2), data = dataset),
      times = B
    )
  )
}

# ---- rmst_fast vs survRM2::rmst2 --------------------------------------------
if (have_survival && have_survRM2) {
  results[["rmst"]] <- summarize_mb(
    "rmst_fast()", "survRM2::rmst2()",
    microbenchmark::microbenchmark(
      fast = rmst_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
                       presorted = TRUE),
      ref  = survRM2::rmst2(time = dataset$tte, status = dataset$event,
                            arm = arm, tau = tau),
      times = B
    )
  )
}

# ---- survdiff_fast(weight = "fh") vs nph::logrank.test ----------------------
if (have_survival && have_nph) {
  results[["wlr"]] <- summarize_mb(
    'survdiff_fast(weight = "fh")', "nph::logrank.test()",
    microbenchmark::microbenchmark(
      fast = survdiff_fast(t_s, e_s, g_s, control = 1, side = 1,
                           weight = "fh", rho = 0, gamma = 1,
                           presorted = TRUE),
      ref  = nph::logrank.test(dataset$tte, dataset$event, dataset$group,
                               rho = 0, gamma = 1),
      times = B
    )
  )
}

# ---- wmst_fast vs survWMST::wmst --------------------------------------------
# survWMST is distributed on GitHub (pauknemj/survWMST), not CRAN.
if (have_survival && have_survWMST) {
  results[["wmst"]] <- summarize_mb(
    "wmst_fast()", "survWMST::wmst()",
    microbenchmark::microbenchmark(
      fast = wmst_fast(t_s, e_s, g_s, control = 1, tau1 = 0, tau2 = tau,
                       side = 1, presorted = TRUE),
      ref  = survWMST::wmst(time = dataset$tte, status = dataset$event,
                            arm = arm, tau0 = 0, tau1 = tau),
      times = B
    )
  )
} else {
  log_line("Skipped wmst_fast benchmark (survWMST not installed).")
}

# ---- milestone_fast vs survfit + summary ------------------------------------
if (have_survival) {
  results[["milestone"]] <- summarize_mb(
    "milestone_fast()", "survfit() + summary()",
    microbenchmark::microbenchmark(
      fast = milestone_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
                            presorted = TRUE),
      ref  = summary(survfit(Surv(tte, event) ~ group, data = dataset),
                     times = tau),
      times = B
    )
  )
}

# ---- medsurv_fast vs nph::nphparams -----------------------------------------
if (have_survival && have_nph) {
  results[["medsurv"]] <- summarize_mb(
    "medsurv_fast()", "nph::nphparams()",
    microbenchmark::microbenchmark(
      fast = medsurv_fast(t_s, e_s, g_s, control = 1, side = 1,
                          method = "nph", presorted = TRUE),
      ref  = nph::nphparams(time = dataset$tte, event = dataset$event,
                            group = as.integer(dataset$group == 2),
                            param_type = "Q", param_par = 0.5),
      times = B
    )
  )
}

# ---- maxcombo_fast vs nph::logrank.maxtest ----------------------------------
if (have_survival && have_nph) {
  results[["maxcombo"]] <- summarize_mb(
    "maxcombo_fast()", "nph::logrank.maxtest()",
    microbenchmark::microbenchmark(
      fast = maxcombo_fast(t_s, e_s, g_s, control = 1, side = 1,
                           rho = c(0, 0, 1), gamma = c(0, 1, 0),
                           presorted = TRUE),
      ref  = nph::logrank.maxtest(dataset$tte, dataset$event,
                                  as.integer(dataset$group == 2)),
      times = B
    )
  )
}

# ---- rmw_fast vs nphRCT::wlrt -----------------------------------------------
# rmw_fast combines a standard and a modestly-weighted log-rank statistic, so
# the reference computes both weighted log-rank components with nphRCT.
if (have_survival && have_nphRCT) {
  results[["rmw"]] <- summarize_mb(
    "rmw_fast()", "nphRCT::wlrt() (two components)",
    microbenchmark::microbenchmark(
      fast = rmw_fast(t_s, e_s, g_s, control = 1, side = 1, s_star = 0.5,
                      presorted = TRUE),
      ref  = {
        nphRCT::wlrt(Surv(tte, event) ~ arm, data = df_rmw,
                     method = "mw", s_star = 1)
        nphRCT::wlrt(Surv(tte, event) ~ arm, data = df_rmw,
                     method = "mw", s_star = 0.5)
      },
      times = B
    )
  )
} else {
  log_line("Skipped rmw_fast benchmark (nphRCT not installed).")
}

# ---- wkm_fast vs nphsim::wkm.Stat -------------------------------------------
# nphsim is distributed on GitHub (keaven/nphsim), not CRAN.
if (have_survival && have_nphsim) {
  results[["wkm"]] <- summarize_mb(
    "wkm_fast()", "nphsim::wkm.Stat()",
    microbenchmark::microbenchmark(
      fast = wkm_fast(t_s, e_s, g_s, control = 1, side = 1, weight = "PF",
                      presorted = TRUE),
      ref  = nphsim::wkm.Stat(survival = dataset$tte,
                              cnsr = 1 - dataset$event,
                              trt = ifelse(dataset$group == 1,
                                           "control", "experimental")),
      times = B
    )
  )
} else {
  log_line("Skipped wkm_fast benchmark (nphsim not installed).")
}

# ---- ahsw_fast vs survAH::ah2 -----------------------------------------------
if (have_survival && have_survAH) {
  results[["ahsw"]] <- summarize_mb(
    "ahsw_fast()", "survAH::ah2()",
    microbenchmark::microbenchmark(
      fast = ahsw_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
                       presorted = TRUE),
      ref  = survAH::ah2(time = dataset$tte, status = dataset$event,
                         arm = arm, tau = tau),
      times = B
    )
  )
} else {
  log_line("Skipped ahsw_fast benchmark (survAH not installed).")
}

# ---- ahr_fast vs AHR::ahrKM --------------------------------------------------
# AHR is archived on CRAN; install with remotes::install_github("cran/AHR").
if (have_survival && have_AHR) {
  results[["ahr"]] <- summarize_mb(
    "ahr_fast()", "AHR::ahrKM()",
    microbenchmark::microbenchmark(
      fast = ahr_fast(t_s, e_s, g_s, control = 1, tau = tau, side = 1,
                      presorted = TRUE),
      ref  = AHR::ahrKM(tau, Surv(tte, event) ~ group, dataset),
      times = B
    )
  )
} else {
  log_line("Skipped ahr_fast benchmark (AHR not installed).")
}

# ---- Collect, print and save -------------------------------------------------
speed_results <- do.call(rbind, results)
rownames(speed_results) <- NULL

cat("\nMedian timings (milliseconds) and speed-up over the reference:\n\n")
print(speed_results, row.names = FALSE)

write.csv(speed_results, file.path(out_dir, "speed_results.csv"),
          row.names = FALSE)
saveRDS(speed_results, file.path(out_dir, "speed_results.rds"))
log_line("Wrote speed_results.csv and speed_results.rds")

cat("\n")
print(sessionInfo())
