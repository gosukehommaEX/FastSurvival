# Compare medsurv_fast() against nph::nphparams() for the non-parametric median
# survival time and the between-group difference, on numerical agreement and on
# speed. Both variance methods of medsurv_fast() are compared: "km" (the native
# Greenwood-plus-kernel-hazard variance) and "nph" (the nph-compatible local
# constant hazard variance). This is a reproducibility script: it uses relative
# paths only, writes a results CSV and a log file, and records the session at
# the end.
#
# Expected behaviour:
#   - Point estimates (control, treatment and difference of medians) agree for
#     both methods, since medsurv_fast() always uses the Kaplan-Meier median.
#   - method = "nph" reproduces the nph::nphparams standard error and p-value to
#     numerical precision whenever the Kaplan-Meier and Nelson-Aalen medians
#     coincide (the usual case).
#   - method = "km" gives a close but not identical standard error, because the
#     quantile variance depends on the local hazard estimate.

log_path <- "tools/benchmarks/compare_medsurv_nphparams_log.txt"
csv_path <- "tools/benchmarks/compare_medsurv_nphparams_results.csv"

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
have_nph <- requireNamespace("nph", quietly = TRUE)
have_microbenchmark <- requireNamespace("microbenchmark", quietly = TRUE)

if (!have_fastsurvival) {
  stop("FastSurvival must be installed (or loaded via devtools::load_all()).")
}
if (!have_nph) {
  log_msg("Package nph is not installed; the comparison cannot run.")
  log_msg("Install it with install.packages(\"nph\") and rerun.")
  quit(save = "no", status = 0)
}

library(FastSurvival)

# Two-piece piecewise-exponential event times: hazard h1 up to the change point
# cp, then hazard h2 afterwards. Sampling is by inversion of the cumulative
# hazard.
rpwexp2 <- function(n, h1, h2, cp) {
  target <- rexp(n, rate = 1)
  threshold <- h1 * cp
  ifelse(target <= threshold, target / h1, cp + (target - threshold) / h2)
}

# Generate a two-group data set with uniform accrual over [0, accrual], an
# administrative censoring time and exponential random dropout. Group 0 is the
# control group, group 1 the treatment group.
gen_data <- function(n_per, rate_c, rate_t = NULL, cp = NULL, hr2 = NULL,
                     accrual = 12, end_time = 36, drop_rate = 0.01, seed = 1) {
  set.seed(seed)
  n_total <- 2 * n_per
  grp <- rep(0:1, each = n_per)
  event_time <- numeric(n_total)
  event_time[grp == 0] <- rexp(n_per, rate = rate_c)
  if (is.null(cp)) {
    event_time[grp == 1] <- rexp(n_per, rate = rate_t)
  } else {
    event_time[grp == 1] <- rpwexp2(n_per, rate_c, rate_c * hr2, cp)
  }
  dropout <- rexp(n_total, rate = drop_rate)
  arrival <- runif(n_total, 0, accrual)
  admin <- end_time - arrival
  cens_time <- pmin(dropout, admin)
  obs_time <- pmin(event_time, cens_time)
  obs_event <- as.integer(event_time <= cens_time)
  data.frame(time = obs_time, event = obs_event, group = grp)
}

# Run both methods on one data set and return a tidy data frame of the
# comparison against nph::nphparams for the median difference and its inference.
compare_one <- function(label, d) {
  ord <- order(d$time)
  np <- nph::nphparams(time = d$time, event = d$event, group = d$group,
                       param_type = "Q", param_par = 0.5)
  nph_val <- c(
    median_control = as.numeric(np$est0),
    median_treatment = as.numeric(np$est1),
    diff = as.numeric(np$tab$Estimate),
    se_diff = as.numeric(np$tab$SE),
    p = as.numeric(np$tab$p_unadj)
  )

  one_method <- function(meth) {
    fast <- medsurv_fast(d$time[ord], d$event[ord], group = d$group[ord],
                         control = 0, side = 2, method = meth, presorted = TRUE)
    fast_val <- c(
      median_control = unname(fast["median.control"]),
      median_treatment = unname(fast["median.treatment"]),
      diff = unname(fast["diff"]),
      se_diff = unname(fast["se.diff"]),
      p = unname(fast["p"])
    )
    data.frame(
      scenario = label,
      method = meth,
      quantity = names(fast_val),
      fast = as.numeric(fast_val),
      nph = as.numeric(nph_val),
      abs_diff = as.numeric(fast_val) - as.numeric(nph_val),
      rel_diff = (as.numeric(fast_val) - as.numeric(nph_val)) /
        ifelse(as.numeric(nph_val) == 0, NA_real_, as.numeric(nph_val)),
      stringsAsFactors = FALSE
    )
  }

  rbind(one_method("km"), one_method("nph"))
}

rate_med <- function(median_time) log(2) / median_time

scenarios <- list(
  list(label = "PH exponential",
       d = gen_data(n_per = 300, rate_c = rate_med(12), rate_t = rate_med(16),
                    seed = 101)),
  list(label = "Delayed effect",
       d = gen_data(n_per = 300, rate_c = rate_med(12), cp = 4, hr2 = 0.55,
                    seed = 202)),
  list(label = "Null (identical)",
       d = gen_data(n_per = 300, rate_c = rate_med(12), rate_t = rate_med(12),
                    seed = 303))
)

results <- do.call(rbind, lapply(scenarios, function(s) compare_one(s$label, s$d)))

log_msg("")
log_msg("Numerical comparison (medsurv_fast vs nph::nphparams, param_type = Q, par = 0.5):")
print(format(results, digits = 4))
write.csv(results, csv_path, row.names = FALSE)

# Summarise the agreement on the standard error of the difference, by method.
se_rows <- results[results$quantity == "se_diff", ]
diff_rows <- results[results$quantity == "diff", ]
log_msg("")
log_msg(paste0("Max abs. difference in the median difference (both methods): ",
               formatC(max(abs(diff_rows$abs_diff)), digits = 4, format = "g")))
for (meth in c("km", "nph")) {
  sub <- se_rows[se_rows$method == meth, ]
  log_msg(paste0("Max rel. difference in SE(diff), method = ", meth, ": ",
                 formatC(max(abs(sub$rel_diff), na.rm = TRUE),
                         digits = 3, format = "g")))
}

# Speed comparison on a larger data set, for both methods against nph.
if (have_microbenchmark) {
  bench_data <- gen_data(n_per = 600, rate_c = rate_med(12), cp = 4, hr2 = 0.55,
                         seed = 404)
  tvec <- bench_data$time
  evec <- bench_data$event
  gvec <- bench_data$group

  # Sort once outside the timing loop, the intended pre-sorted fast path.
  ord_b <- order(tvec)
  tvs <- tvec[ord_b]; evs <- evec[ord_b]; gvs <- gvec[ord_b]

  timing <- microbenchmark::microbenchmark(
    fast_km = medsurv_fast(tvs, evs, group = gvs, control = 0, side = 2,
                           method = "km", presorted = TRUE),
    fast_nph = medsurv_fast(tvs, evs, group = gvs, control = 0, side = 2,
                            method = "nph", presorted = TRUE),
    nph = nph::nphparams(time = tvec, event = evec, group = gvec,
                         param_type = "Q", param_par = 0.5),
    times = 100
  )
  log_msg("")
  log_msg("Speed comparison (n = 600 per group, 100 replicates):")
  print(summary(timing))

  med_km <- median(timing$time[timing$expr == "fast_km"])
  med_nph_fast <- median(timing$time[timing$expr == "fast_nph"])
  med_nph <- median(timing$time[timing$expr == "nph"])
  log_msg(paste0("Median speedup factor (nph / fast_km): ",
                 formatC(med_nph / med_km, digits = 3, format = "g")))
  log_msg(paste0("Median speedup factor (nph / fast_nph): ",
                 formatC(med_nph / med_nph_fast, digits = 3, format = "g")))
} else {
  log_msg("Package microbenchmark is not installed; the speed comparison was skipped.")
}

log_msg("")
log_msg("Session information:")
si <- capture.output(sessionInfo())
cat(si, sep = "\n", file = log_path, append = TRUE)
cat(si, sep = "\n")
cat("\n")
