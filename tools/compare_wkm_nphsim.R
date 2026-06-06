# Compare wkm_fast() against the nphsim weighted Kaplan-Meier (Pepe-Fleming)
# test, on numerical agreement and on speed.
#
# The nphsim package (keaven/nphsim, GitHub only) depends on survMisc, which was
# archived from CRAN on 2026-03-17, so a plain install of nphsim currently
# fails. To keep the comparison runnable, this script bundles a faithful
# reproduction of nphsim::wkm.Stat that needs only the survival package, and
# uses it when nphsim is not installed. If nphsim is available, the real
# nphsim::wkm.Stat is used instead. Either way wkm_fast(side = 1, weight = "PF")
# is expected to reproduce the reference z statistic and one-sided p-value to
# numerical precision for data without tied event times, and to be much faster.
#
# This is a reproducibility script: relative paths only, a results CSV and a log
# file, and the session is recorded at the end.

log_path <- "tools/benchmarks/compare_wkm_nphsim_log.txt"
csv_path <- "tools/benchmarks/compare_wkm_nphsim_results.csv"

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
have_survival <- requireNamespace("survival", quietly = TRUE)
have_nphsim <- requireNamespace("nphsim", quietly = TRUE)
have_microbenchmark <- requireNamespace("microbenchmark", quietly = TRUE)

if (!have_fastsurvival) {
  stop("FastSurvival must be installed (or loaded via devtools::load_all()).")
}
if (!have_survival) {
  stop("The survival package is required for the reference computation.")
}

library(FastSurvival)

# Faithful reproduction of nphsim::wkm.Stat (keaven/nphsim), the Pepe-Fleming
# (1989) weighted Kaplan-Meier test, using only the survival package. Provided
# so the comparison runs without installing nphsim. group_treat: 1 = treatment
# (experimental), 0 = control. Returns the one-sided z and p-value for treatment
# superiority, matching the orientation of nphsim::wkm.Stat.
wkm_stat_ref <- function(time, event, group_treat) {
  n <- length(time)
  n1 <- sum(group_treat == 1)
  n2 <- sum(group_treat == 0)
  d1 <- data.frame(t = time[group_treat == 1], s = event[group_treat == 1])
  d2 <- data.frame(t = time[group_treat == 0], s = event[group_treat == 0])

  f1 <- summary(survival::survfit(survival::Surv(t, s == 1) ~ 1, data = d1))
  f2 <- summary(survival::survfit(survival::Surv(t, s == 1) ~ 1, data = d2))
  s1_fun <- approxfun(f1$time, f1$surv, method = "constant",
                      yleft = 1, rule = 2, f = 0)
  s2_fun <- approxfun(f2$time, f2$surv, method = "constant",
                      yleft = 1, rule = 2, f = 0)

  et <- sort(time)
  grid <- (c(0, et[-n]) + et) / 2
  width <- diff(c(0, et))
  km_s1 <- s1_fun(grid)
  km_s2 <- s2_fun(grid)

  c1 <- summary(survival::survfit(survival::Surv(t, s == 0) ~ 1, data = d1))
  c2 <- summary(survival::survfit(survival::Surv(t, s == 0) ~ 1, data = d2))
  g1_fun <- if (length(c1$time) >= 1) {
    approxfun(c(c1$time, max(c1$time) + 1), c(1, c1$surv), method = "constant",
              yleft = 1, yright = min(c1$surv), rule = 2, f = 1)
  } else {
    function(x) rep(1, length(x))
  }
  g2_fun <- if (length(c2$time) >= 1) {
    approxfun(c(c2$time, max(c2$time) + 1), c(1, c2$surv), method = "constant",
              yleft = 1, yright = min(c2$surv), rule = 2, f = 1)
  } else {
    function(x) rep(1, length(x))
  }
  g1 <- g1_fun(grid)
  g2 <- g2_fun(grid)
  wt <- ifelse(g1 + g2 == 0, 0, (n * g1 * g2) / (n1 * g1 + n2 * g2))

  num_raw <- sum(wt * (km_s1 - km_s2) * width)

  fp <- summary(survival::survfit(survival::Surv(time, event) ~ 1))
  sp_fun <- approxfun(fp$time, fp$surv, method = "constant",
                      yleft = 1, rule = 2, f = 0)
  sp <- sp_fun(grid)
  a_seq <- cumsum(width * wt * sp)
  a_rem <- a_seq[n] - a_seq
  sp_next <- c(sp[-1], min(fp$surv))
  d_sm <- sp_next - sp
  keep <- wt > 0 & sp > 0
  variance <- -sum((a_rem[keep]^2 / (sp[keep]^2 * wt[keep])) * d_sm[keep])

  z <- sqrt(n1 * n2 / n) * num_raw / sqrt(variance)
  list(z = z, pval = 1 - pnorm(z))
}

# Choose the reference: the real nphsim if installed, otherwise the bundled
# survival-only reproduction.
if (have_nphsim) {
  ref_label <- "nphsim::wkm.Stat"
  ref_fun <- function(d) {
    trt <- ifelse(d$group == 1, "experimental", "control")
    np <- nphsim::wkm.Stat(survival = d$time, cnsr = 1 - d$event, trt = trt)
    list(z = as.numeric(np$z), pval = as.numeric(np$pval))
  }
} else {
  ref_label <- "survival-based reproduction of nphsim::wkm.Stat"
  ref_fun <- function(d) wkm_stat_ref(d$time, d$event, d$group)
}
log_msg(paste("Reference implementation:", ref_label))

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

compare_one <- function(label, d) {
  fast <- wkm_fast(d$time, d$event, group = d$group, control = 0, side = 1,
                   weight = "PF")
  ref <- ref_fun(d)

  fast_val <- c(z = unname(fast["z"]), p = unname(fast["p"]))
  ref_val <- c(z = as.numeric(ref$z), p = as.numeric(ref$pval))

  data.frame(
    scenario = label,
    quantity = names(fast_val),
    fast = as.numeric(fast_val),
    reference = as.numeric(ref_val),
    abs_diff = as.numeric(fast_val) - as.numeric(ref_val),
    rel_diff = (as.numeric(fast_val) - as.numeric(ref_val)) /
      ifelse(as.numeric(ref_val) == 0, NA_real_, as.numeric(ref_val)),
    stringsAsFactors = FALSE
  )
}

rate_med <- function(median_time) log(2) / median_time

scenarios <- list(
  list(label = "PH exponential",
       d = gen_data(n_per = 300, rate_c = rate_med(12), rate_t = rate_med(16),
                    seed = 101)),
  list(label = "Delayed effect",
       d = gen_data(n_per = 300, rate_c = rate_med(12), cp = 4, hr1 = 1,
                    hr2 = 0.55, seed = 202)),
  list(label = "Crossing hazards",
       d = gen_data(n_per = 300, rate_c = rate_med(12), cp = 6, hr1 = 1.6,
                    hr2 = 0.5, seed = 303)),
  list(label = "Null (identical)",
       d = gen_data(n_per = 300, rate_c = rate_med(12), rate_t = rate_med(12),
                    seed = 404))
)

results <- do.call(rbind, lapply(scenarios, function(s) compare_one(s$label, s$d)))

log_msg("")
log_msg("Numerical comparison (wkm_fast side = 1, weight = PF vs reference):")
print(format(results, digits = 6))
write.csv(results, csv_path, row.names = FALSE)

z_rows <- results[results$quantity == "z", ]
p_rows <- results[results$quantity == "p", ]
log_msg("")
log_msg(paste0("Max abs. difference in z: ",
               formatC(max(abs(z_rows$abs_diff)), digits = 4, format = "g")))
log_msg(paste0("Max abs. difference in one-sided p: ",
               formatC(max(abs(p_rows$abs_diff)), digits = 4, format = "g")))

if (have_microbenchmark) {
  bench_data <- gen_data(n_per = 600, rate_c = rate_med(12), cp = 4, hr1 = 1,
                         hr2 = 0.55, seed = 505)
  tvec <- bench_data$time
  evec <- bench_data$event
  gvec <- bench_data$group

  timing <- microbenchmark::microbenchmark(
    fast = wkm_fast(tvec, evec, group = gvec, control = 0, side = 1,
                    weight = "PF"),
    reference = ref_fun(bench_data),
    times = 100
  )
  log_msg("")
  log_msg("Speed comparison (n = 600 per group, 100 replicates):")
  print(summary(timing))

  med_fast <- median(timing$time[timing$expr == "fast"])
  med_ref <- median(timing$time[timing$expr == "reference"])
  log_msg(paste0("Median speedup factor (reference / fast): ",
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
