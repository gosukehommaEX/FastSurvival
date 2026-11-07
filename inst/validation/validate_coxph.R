# Validation of coxph_fast() and of the "coxph" statistic of analysis_fast()
# against survival::coxph().
#
# Run from the package root with the version to be validated installed (for
# example after devtools::install()) or loaded with devtools::load_all():
#   source("inst/validation/validate_coxph.R")
# The run takes a few minutes. It writes to inst/validation/results/:
#   validate_coxph_summary.md   summary, generated from the results below
#   validate_coxph_summary.csv  numbers of fits by group and verdict
#   validate_coxph_nonpass.csv  every fit that is not PASS
#
# Every fit is compared with survival::coxph() with the same handling of tied
# event times and coxph.control(timefix = FALSE), so that, as in coxph_fast(),
# only exactly equal times are tied. Each fit receives one verdict:
#   PASS       both converge, the log hazard ratios and the standard errors
#              differ by at most 1e-10 relative to max(1, |coxph() value|),
#              and, for coxph_fast(), the numbers of iterations are equal.
#   EXPLAINED  one of the registered reasons below applies.
#   FAIL       otherwise.
# Registered reasons:
#   E1  The estimate is infinite: coxph_fast() returns NA, and coxph() gives
#       a warning, a coefficient that is not finite, or a standard error
#       above 1e3.
#   E2  Halved steps: both converge, the log hazard ratios differ by at most
#       1e-6, the iterations differ (coxph_fast() keeps a bracket of the root
#       after a halved step, see ?coxph_fast), and the score test of coxph()
#       at the coxph_fast() estimate is at most 1e-8.
#   E3  coxph() does not converge in 20 iterations while coxph_fast() does,
#       and the score test of coxph() at the coxph_fast() estimate is at most
#       1e-8.
#   E4  coxph() warns that the coefficient may be infinite, but both
#       converge, the log hazard ratios differ by at most 1e-6, and the score
#       test of coxph() at the coxph_fast() estimate is at most 1e-8 (the
#       estimate is finite, so the warning is a false alarm).
# A self-test changes the coxph_fast() estimates of the PASS fits by 1e-8
# times max(1, |estimate|) and checks that none of them is then classified
# as PASS.

library(survival)
if (!"FastSurvival" %in% loadedNamespaces()) library(FastSurvival)
stopifnot("ties" %in% names(formals(FastSurvival::coxph_fast)))

res_dir <- file.path("inst", "validation", "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

tol_pass  <- 1e-10
tol_e2    <- 1e-6
tol_score <- 1e-8

# Reference fit by survival::coxph(); warnings are kept as text. With
# iter.max = 0 the fit stays at init, and its score test is at init.
ref_fit <- function(time, event, x, strata, ties, init = NULL,
                    iter.max = 20) {
  d <- data.frame(time = time, event = event, x = x)
  f <- Surv(time, event) ~ x
  if (!is.null(strata)) {
    d$s <- strata
    f <- Surv(time, event) ~ x + strata(s)
  }
  ctrl <- coxph.control(timefix = FALSE, iter.max = iter.max)
  warn <- character(0)
  fit <- withCallingHandlers(
    if (is.null(init)) {
      coxph(f, data = d, ties = ties, control = ctrl)
    } else {
      coxph(f, data = d, ties = ties, init = init, control = ctrl)
    },
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(coef = unname(fit$coefficients[1L]), se = sqrt(fit$var[1L, 1L]),
       iter = fit$iter,
       score = if (is.null(fit$score)) NA_real_ else unname(fit$score),
       warning = paste(warn, collapse = " | "))
}

rel_diff <- function(a, b) abs(a - b) / pmax(1, abs(b))

# Verdicts of a data frame of fits; iter_fast is NA for analysis_fast()
classify <- function(d) {
  conv <- d$converged
  ok_ref <- d$warning_cox == "" & is.finite(d$coef_cox)
  same_iter <- is.na(d$iter_fast) | d$iter_fast == d$iter_cox
  pass <- !is.na(conv) & conv & ok_ref &
    rel_diff(d$coef_fast, d$coef_cox) <= tol_pass &
    rel_diff(d$se_fast, d$se_cox) <= tol_pass & same_iter
  pass[is.na(pass)] <- FALSE
  score_ok <- !is.na(d$score_at_fast) & d$score_at_fast <= tol_score
  e1 <- is.na(conv) & (d$warning_cox != "" | !is.finite(d$coef_cox) |
                         (!is.na(d$se_cox) & d$se_cox > 1e3))
  e2 <- !pass & !is.na(conv) & conv & ok_ref & score_ok &
    abs(d$coef_fast - d$coef_cox) <= tol_e2 & !same_iter
  e3 <- !pass & !is.na(conv) & conv & score_ok &
    grepl("Ran out of iterations", d$warning_cox)
  e4 <- !pass & !is.na(conv) & conv & score_ok &
    grepl("may be infinite", d$warning_cox) &
    !grepl("Ran out of iterations", d$warning_cox) &
    abs(d$coef_fast - d$coef_cox) <= tol_e2
  e1[is.na(e1)] <- FALSE
  e2[is.na(e2)] <- FALSE
  e3[is.na(e3)] <- FALSE
  e4[is.na(e4)] <- FALSE
  reason <- ifelse(pass, "", ifelse(e1, "E1", ifelse(e2, "E2",
                   ifelse(e3, "E3", ifelse(e4, "E4", "")))))
  verdict <- ifelse(pass, "PASS", ifelse(reason != "", "EXPLAINED", "FAIL"))
  data.frame(verdict = verdict, reason = reason)
}

# The score test of coxph() at the estimate, computed only when needed
score_needed <- function(conv, cf, sf, itf, ref) {
  isTRUE(conv) && !(ref$warning == "" && is.finite(ref$coef) &&
                      rel_diff(cf, ref$coef) <= tol_pass &&
                      rel_diff(sf, ref$se) <= tol_pass &&
                      (is.na(itf) || itf == ref$iter))
}

one_fit <- function(time, event, x, strata, ties) {
  fast <- suppressWarnings(
    coxph_fast(time, event, x, control = 0, strata = strata, ties = ties))
  ref <- ref_fit(time, event, x, strata, ties)
  conv <- attr(fast, "converged")
  cf <- unname(fast["coef"])
  sf <- unname(fast["se(coef)"])
  itf <- attr(fast, "iter")
  score <- NA_real_
  if (score_needed(conv, cf, sf, itf, ref)) {
    score <- ref_fit(time, event, x, strata, ties, init = cf,
                     iter.max = 0)$score
  }
  data.frame(coef_fast = cf, se_fast = sf, iter_fast = itf, converged = conv,
             coef_cox = ref$coef, se_cox = ref$se, iter_cox = ref$iter,
             warning_cox = ref$warning, score_at_fast = score)
}

# Two-group data: exponential event and censoring times, optionally rounded
# to quarters (tied times) and with three strata of different baseline hazards
sim_one <- function(n, hr, p_trt, cens_rate, round_time, n_strata) {
  x <- as.integer(runif(n) < p_trt)
  s <- if (n_strata > 1L) sample.int(n_strata, n, replace = TRUE) else NULL
  base <- if (is.null(s)) 1 else c(1, 0.5, 2)[s]
  ev <- rexp(n, base * ifelse(x == 1L, hr, 1))
  ce <- rexp(n, cens_rate)
  time <- pmin(ev, ce)
  if (round_time) time <- ceiling(time * 4) / 4
  list(time = time, event = as.integer(ev <= ce), x = x, strata = s)
}

run_grid <- function(group, grid, n_rep) {
  rows <- list()
  for (g in seq_len(nrow(grid))) {
    p <- grid[g, ]
    for (rep in seq_len(n_rep)) {
      d <- sim_one(p$n, p$hr, p$p_trt, p$cens_rate, p$round_time,
                   p$n_strata)
      if (length(unique(d$x)) < 2L || sum(d$event) == 0L) next
      for (ties in if (p$round_time) c("efron", "breslow") else "efron") {
        rows[[length(rows) + 1L]] <- data.frame(
          group = group, scenario = g, rep = rep, as.list(p), ties = ties,
          events = sum(d$event), n_trt = sum(d$x),
          one_fit(d$time, d$event, d$x, d$strata, ties))
      }
    }
  }
  do.call(rbind, rows)
}

set.seed(20261108)
message("Group A (regular data) ...")
grid_a <- expand.grid(n = c(30L, 100L, 500L, 2000L),
                      hr = c(0.5, 0.8, 1, 1.25, 2), p_trt = c(0.5, 0.25),
                      cens_rate = c(0.2, 1), round_time = c(FALSE, TRUE),
                      n_strata = c(1L, 3L), KEEP.OUT.ATTRS = FALSE)
res_a <- run_grid("A regular", grid_a, n_rep = 15L)
message("Group B (small treatment group, extreme hazard ratio) ...")
grid_b <- expand.grid(n = c(20L, 100L), hr = c(0.01, 100), p_trt = 0.05,
                      cens_rate = c(0.2, 2), round_time = c(FALSE, TRUE),
                      n_strata = c(1L, 3L), KEEP.OUT.ATTRS = FALSE)
res_b <- run_grid("B extreme", grid_b, n_rep = 25L)

# ---- Group C: the "coxph" statistic of analysis_fast() --------------------
# Each row of analysis_fast() is compared with coxph() on the data cut at the
# calendar time of the row (subjects enrolled by then; events after it are
# censored at it).
cut_data <- function(dat, sim_id, cutoff) {
  d <- dat[dat$sim == sim_id & dat$accrual_time <= cutoff, ]
  before <- (d$accrual_time + d$tte) <= cutoff
  t <- ifelse(before, d$tte, cutoff - d$accrual_time)
  e <- ifelse(before, d$event, 0L)
  list(time = t, event = as.integer(e), x = as.integer(d$group != 1),
       subgroup = if ("subgroup" %in% names(d)) d$subgroup)
}

message("Group C (analysis_fast) ...")
dat_c1 <- simdata_fast(nsim = 100, n = c(150, 150), a.time = c(0, 12),
                       a.prop = 1,
                       e.hazard = list(log(2) / 12, log(2) / 16),
                       d.median = list(36, 36), seed = 11)
dat_c2 <- dat_c1
dat_c2$tte <- ceiling(dat_c2$tte)
dat_c3 <- simdata_fast(nsim = 100, n = c(150, 150), a.time = c(0, 12),
                       a.prop = 1,
                       e.hazard = list(list(0.10, 0.05), list(0.07, 0.035)),
                       prevalence = c(0.5, 0.5), seed = 13)
configs_c <- list(
  list(label = "continuous times, calendar looks", dat = dat_c1,
       ties = "efron", strata = NULL, looks = list(time.looks = c(12, 24))),
  list(label = "continuous times, event looks", dat = dat_c1,
       ties = "efron", strata = NULL, looks = list(event.looks = c(80, 160))),
  list(label = "whole-month times, Efron", dat = dat_c2, ties = "efron",
       strata = NULL, looks = list(time.looks = c(12, 24))),
  list(label = "whole-month times, Breslow", dat = dat_c2, ties = "breslow",
       strata = NULL, looks = list(time.looks = c(12, 24))),
  list(label = "stratified by subgroup", dat = dat_c3, ties = "efron",
       strata = "subgroup", looks = list(time.looks = c(12, 24)))
)
rows <- list()
for (k in seq_along(configs_c)) {
  cfg <- configs_c[[k]]
  res <- do.call(analysis_fast,
                 c(list(cfg$dat, control = 1, stat = "coxph",
                        ties = cfg$ties, strata = cfg$strata), cfg$looks))
  for (r in seq_len(nrow(res))) {
    if (!isTRUE(res$reached[r])) next
    cd <- cut_data(cfg$dat, res$sim[r], res$cutoff[r])
    if (length(unique(cd$x)) < 2L || sum(cd$event) == 0L) next
    st <- if (is.null(cfg$strata)) NULL else cd$subgroup
    ref <- ref_fit(cd$time, cd$event, cd$x, st, cfg$ties)
    cf <- res$cox.coef[r]
    sf <- res$cox.se[r]
    conv <- if (is.na(cf)) NA else TRUE
    score <- NA_real_
    if (score_needed(conv, cf, sf, NA, ref)) {
      score <- ref_fit(cd$time, cd$event, cd$x, st, cfg$ties, init = cf,
                       iter.max = 0)$score
    }
    rows[[length(rows) + 1L]] <- data.frame(
      group = "C analysis_fast", scenario = k, rep = res$sim[r],
      config = cfg$label, look = res$look[r], ties = cfg$ties,
      events = sum(cd$event), n_trt = sum(cd$x),
      coef_fast = cf, se_fast = sf, iter_fast = NA_integer_,
      converged = conv, coef_cox = ref$coef, se_cox = ref$se,
      iter_cox = ref$iter, warning_cox = ref$warning, score_at_fast = score)
  }
}
res_c <- do.call(rbind, rows)

# ---- Verdicts, self-test, and output ---------------------------------------
common <- c("group", "scenario", "rep", "ties", "events", "n_trt",
            "coef_fast", "se_fast", "iter_fast", "converged", "coef_cox",
            "se_cox", "iter_cox", "warning_cox", "score_at_fast")
all_fits <- rbind(res_a[common], res_b[common], res_c[common])
all_fits <- cbind(all_fits, classify(all_fits))

pass <- all_fits$verdict == "PASS"
perturbed <- all_fits[pass, ]
perturbed$coef_fast <- perturbed$coef_fast +
  1e-8 * pmax(1, abs(perturbed$coef_fast))
self_v <- classify(perturbed)$verdict
self_ok <- sum(pass) > 0L && !any(self_v == "PASS")

groups <- unique(all_fits$group)
summary_tab <- do.call(rbind, lapply(groups, function(g) {
  v <- all_fits$verdict[all_fits$group == g]
  data.frame(group = g, fits = length(v), PASS = sum(v == "PASS"),
             EXPLAINED = sum(v == "EXPLAINED"), FAIL = sum(v == "FAIL"))
}))
write.csv(summary_tab, file.path(res_dir, "validate_coxph_summary.csv"),
          row.names = FALSE)
write.csv(all_fits[!pass, ], file.path(res_dir, "validate_coxph_nonpass.csv"),
          row.names = FALSE)

fmt <- function(x) format(signif(x, 2))
reasons <- table(factor(all_fits$reason[all_fits$verdict == "EXPLAINED"],
                        levels = c("E1", "E2", "E3", "E4")))
md <- c(
  "# Validation of coxph_fast() against survival::coxph()",
  "",
  paste0("Generated by `inst/validation/validate_coxph.R` on ",
         format(Sys.Date()), " with ", R.version.string, ", FastSurvival ",
         as.character(utils::packageVersion("FastSurvival")),
         ", and survival ", as.character(utils::packageVersion("survival")),
         ". Do not edit by hand."),
  "",
  "| Group | Fits | PASS | EXPLAINED | FAIL |",
  "|---|---|---|---|---|",
  sprintf("| %s | %d | %d | %d | %d |", summary_tab$group, summary_tab$fits,
          summary_tab$PASS, summary_tab$EXPLAINED, summary_tab$FAIL),
  "",
  paste0("Group A: ", nrow(grid_a), " scenarios with 15 data sets each ",
         "(n 30 to 2000, hazard ratio 0.5 to 2, 25% or 50% treated, two ",
         "censoring rates, continuous or tied times, 1 or 3 strata). ",
         "Group B: ", nrow(grid_b), " scenarios with 25 data sets each ",
         "(n 20 or 100, 5% treated, hazard ratio 0.01 or 100). Data sets ",
         "with one group only or no events are skipped. Tied times are ",
         "fitted with both the Efron and the Breslow handling. Group C: ",
         "the rows of analysis_fast() for ", length(configs_c),
         " configurations of 100 simulated trials and two looks."),
  "",
  paste0("Largest differences among the PASS fits (relative to ",
         "max(1, |coxph() value|)): log hazard ratio ",
         fmt(max(rel_diff(all_fits$coef_fast[pass], all_fits$coef_cox[pass]))),
         ", standard error ",
         fmt(max(rel_diff(all_fits$se_fast[pass], all_fits$se_cox[pass]))),
         " (tolerance ", fmt(tol_pass), ")."),
  "",
  paste0("EXPLAINED: E1 (infinite estimate) ", reasons[["E1"]],
         ", E2 (halved steps) ", reasons[["E2"]],
         ", E3 (coxph() did not converge) ", reasons[["E3"]],
         ", E4 (false warning of an infinite coefficient by coxph()) ",
         reasons[["E4"]],
         ". The reasons are defined at the top of the script; every fit ",
         "that is not PASS is listed in `validate_coxph_nonpass.csv`."),
  "",
  paste0("Self-test: after changing the ", sum(pass), " PASS estimates by ",
         "1e-8 times max(1, |estimate|), ", sum(self_v == "PASS"),
         " are PASS and ",
         sum(self_v == "FAIL"), " are FAIL (",
         if (self_ok) "the comparison detects the change" else
           "THE COMPARISON DOES NOT DETECT THE CHANGE", ").")
)
writeLines(md, file.path(res_dir, "validate_coxph_summary.md"))
message(paste(md, collapse = "\n"))
if (any(all_fits$verdict == "FAIL") || !self_ok) {
  warning("validate_coxph: FAIL in ", sum(all_fits$verdict == "FAIL"),
          " fits or the self-test failed; see ", res_dir)
}
