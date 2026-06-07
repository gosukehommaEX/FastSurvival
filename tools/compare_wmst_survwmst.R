# Compare wmst_fast() against survWMST::wmst (pauknemj/survWMST) for the window
# mean survival time, on numerical agreement and on speed. survWMST is a GitHub
# package and is therefore not used in the automated test suite; this script is
# meant to be run manually where survWMST has been installed. It is a
# reproducibility script: relative paths only, a results CSV and a log file,
# and the session is recorded at the end.
#
# Argument mapping: survWMST::wmst uses tau0 for the window start and tau1 for
# the window end, whereas wmst_fast uses tau1 for the window start and tau2 for
# the window end. For a window [a, b] this script calls wmst_fast(tau1 = a,
# tau2 = b) and survWMST::wmst(tau0 = a, tau1 = b).
#
# Expected behaviour: wmst_fast reproduces the per-arm WMST, the difference, its
# standard error and the two-sided p-value of survWMST::wmst to numerical
# precision, and is substantially faster.

log_path <- "tools/benchmarks/compare_wmst_survwmst_log.txt"
csv_path <- "tools/benchmarks/compare_wmst_survwmst_results.csv"

if (!dir.exists("tools/benchmarks")) {
  dir.create("tools/benchmarks", recursive = TRUE, showWarnings = FALSE)
}

log_msg <- function(text) {
  cat(text, "\n", file = log_path, append = TRUE)
  cat(text, "\n")
}

cat("", file = log_path, append = FALSE)
log_msg(paste("Comparison run:", format(Sys.time())))

have_fastsurvival <- requireNamespace("FastSurvival", quietly = TRUE)
have_survwmst <- requireNamespace("survWMST", quietly = TRUE)
have_microbenchmark <- requireNamespace("microbenchmark", quietly = TRUE)

if (!have_fastsurvival) {
  stop("FastSurvival must be installed (or loaded via devtools::load_all()).")
}
if (!have_survwmst) {
  log_msg("Package survWMST is not installed; the comparison cannot run.")
  log_msg("Install it with remotes::install_github(\"pauknemj/survWMST\") and rerun.")
  quit(save = "no", status = 0)
}

library(FastSurvival)

# Two-piece piecewise-exponential event times.
rpwexp2 <- function(n, h1, h2, cp) {
  target <- rexp(n, rate = 1)
  threshold <- h1 * cp
  ifelse(target <= threshold, target / h1, cp + (target - threshold) / h2)
}

# Generate a two-group data set. Group 0 is control, group 1 treatment. The
# treatment group is exponential (rate_t) or two-piece piecewise-exponential
# with early hazard ratio hr1 and late hazard ratio hr2 relative to control.
gen_data <- function(n_per, rate_c, rate_t = NULL, cp = NULL, hr1 = 1, hr2 = NULL,
                     accrual = 12, end_time = 36, drop_rate = 0.01, seed = 1) {
  set.seed(seed)
  n_total <- 2 * n_per
  grp <- rep(0:1, each = n_per)
  event_time <- numeric(n_total)
  event_time[grp == 0] <- rexp(n_per, rate = rate_c)
  if (is.null(cp)) {
    event_time[grp == 1] <- rexp(n_per, rate = rate_t)
  } else {
    event_time[grp == 1] <- rpwexp2(n_per, rate_c * hr1, rate_c * hr2, cp)
  }
  dropout <- rexp(n_total, rate = drop_rate)
  arrival <- runif(n_total, 0, accrual)
  admin <- end_time - arrival
  cens_time <- pmin(dropout, admin)
  obs_time <- pmin(event_time, cens_time)
  obs_event <- as.integer(event_time <= cens_time)
  data.frame(time = obs_time, event = obs_event, group = grp)
}

# Coerce any value to a length-one numeric, returning NA when empty or missing.
safe_num <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0) NA_real_ else v[[1]]
}

# Extract the per-arm WMST estimate and its standard error from the named
# vector ref$WMST.arm{k}$wmst (elements "Est." and "se"). Returns c(est, se).
arm_wmst <- function(ref, k) {
  v <- ref[[paste0("WMST.arm", k)]]$wmst
  c(safe_num(v["Est."]), safe_num(v["se"]))
}

# Run wmst_fast and survWMST::wmst on one data set and window [w_start, w_end],
# and return a tidy data frame comparing the per-arm WMST, the difference, its
# standard error and the two-sided p-value. The per-arm estimates and standard
# errors come from ref$WMST.arm{k}$wmst; the difference and p-value come from
# the between-group contrast table ref$result (row 1, the difference row); the
# standard error of the difference is the root sum of squared per-arm standard
# errors, as survWMST computes it.
compare_one <- function(label, d, w_start, w_end) {
  fast <- wmst_fast(d$time, d$event, group = d$group, control = 0,
                    tau1 = w_start, tau2 = w_end, side = 2)
  ref <- suppressWarnings(survWMST::wmst(time = d$time, status = d$event,
                                         arm = d$group,
                                         tau0 = w_start, tau1 = w_end))

  a1 <- arm_wmst(ref, 1)
  a0 <- arm_wmst(ref, 0)
  arm1_est <- a1[1]
  arm0_est <- a0[1]
  se_diff_ref <- sqrt(a1[2]^2 + a0[2]^2)

  contrast <- ref$result
  diff_ref <- safe_num(contrast[1, 1])
  p_ref <- safe_num(contrast[1, ncol(contrast)])

  fast_val <- c(
    wmst_control = unname(fast["wmst.control"]),
    wmst_treatment = unname(fast["wmst.treatment"]),
    diff = unname(fast["diff"]),
    se_diff = unname(fast["se.diff"]),
    p = unname(fast["p"])
  )
  ref_val <- c(
    wmst_control = arm0_est,
    wmst_treatment = arm1_est,
    diff = diff_ref,
    se_diff = se_diff_ref,
    p = p_ref
  )

  data.frame(
    scenario = label,
    window = sprintf("[%g, %g]", w_start, w_end),
    quantity = names(fast_val),
    fast = as.numeric(fast_val),
    survwmst = as.numeric(ref_val),
    abs_diff = as.numeric(fast_val) - as.numeric(ref_val),
    rel_diff = (as.numeric(fast_val) - as.numeric(ref_val)) /
      ifelse(as.numeric(ref_val) == 0, NA_real_, as.numeric(ref_val)),
    stringsAsFactors = FALSE
  )
}

rate_med <- function(median_time) log(2) / median_time

scenarios <- list(
  list(label = "PH exponential", w_start = 0, w_end = 12,
       d = gen_data(n_per = 300, rate_c = rate_med(12), rate_t = rate_med(16),
                    seed = 101)),
  list(label = "Delayed effect", w_start = 4, w_end = 14,
       d = gen_data(n_per = 300, rate_c = rate_med(12), cp = 4, hr1 = 1,
                    hr2 = 0.55, seed = 202)),
  list(label = "Crossing hazards", w_start = 6, w_end = 14,
       d = gen_data(n_per = 300, rate_c = rate_med(12), cp = 6, hr1 = 1.6,
                    hr2 = 0.5, seed = 303)),
  list(label = "Null (identical)", w_start = 0, w_end = 12,
       d = gen_data(n_per = 300, rate_c = rate_med(12), rate_t = rate_med(12),
                    seed = 404))
)

results <- do.call(rbind, lapply(scenarios, function(s) {
  compare_one(s$label, s$d, s$w_start, s$w_end)
}))

log_msg("")
log_msg("Numerical comparison (wmst_fast vs survWMST::wmst):")
print(format(results, digits = 6))
write.csv(results, csv_path, row.names = FALSE)

log_msg("")
for (q in c("diff", "se_diff", "p")) {
  sub <- results[results$quantity == q, ]
  log_msg(paste0("Max abs. difference in ", q, ": ",
                 formatC(max(abs(sub$abs_diff)), digits = 4, format = "g")))
}

# Speed comparison on a larger data set over a window.
if (have_microbenchmark) {
  bench_data <- gen_data(n_per = 600, rate_c = rate_med(12), cp = 4, hr1 = 1,
                         hr2 = 0.55, seed = 505)
  tvec <- bench_data$time
  evec <- bench_data$event
  gvec <- bench_data$group

  timing <- microbenchmark::microbenchmark(
    fast = wmst_fast(tvec, evec, group = gvec, control = 0,
                     tau1 = 4, tau2 = 14, side = 2),
    survwmst = suppressWarnings(survWMST::wmst(time = tvec, status = evec,
                                               arm = gvec, tau0 = 4, tau1 = 14)),
    times = 100
  )
  log_msg("")
  log_msg("Speed comparison (n = 600 per group, 100 replicates):")
  print(summary(timing))

  med_fast <- median(timing$time[timing$expr == "fast"])
  med_ref <- median(timing$time[timing$expr == "survwmst"])
  log_msg(paste0("Median speedup factor (survWMST / fast): ",
                 formatC(med_ref / med_fast, digits = 3, format = "g")))
} else {
  log_msg("Package microbenchmark is not installed; the speed comparison was skipped.")
}

log_msg("")
log_msg("Session information:")
si <- capture.output(sessionInfo())
cat(si, sep = "\n", file = log_path, append = TRUE)
cat(si, sep = "\n")
cat("\n")
