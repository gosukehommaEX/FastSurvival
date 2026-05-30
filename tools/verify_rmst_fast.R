# Verification of rmst_fast() against survRM2::rmst2() and survival.
# Reproducibility script kept as an artifact (not part of the package).
# During development, load the package first, e.g. devtools::load_all().
# All paths are relative; the log is written next to this script.

stopifnot(exists("rmst_fast"))
suppressMessages(library(survival))

has_survRM2 <- requireNamespace("survRM2", quietly = TRUE)

log_path <- "rmst_fast_verification_log.txt"
cat("rmst_fast verification log\n", file = log_path, append = FALSE)
cat(format(Sys.time()), "\n", file = log_path, append = TRUE)
cat("survRM2 available: ", has_survRM2, "\n\n", file = log_path, append = TRUE)

# Reference 1: survival's restricted mean via survmean (internal but stable).
rmst_survival <- function(tt, ss, tau) {
  sf <- survfit(Surv(tt, ss) ~ 1)
  mm <- tryCatch(survival:::survmean(sf, rmean = tau)$matrix,
                 error = function(e) NULL)
  if (is.null(mm)) return(c(est = NA_real_, se = NA_real_))
  est_nm <- names(mm)[grepl("rmean", names(mm)) & !grepl("se", names(mm))]
  se_nm  <- names(mm)[grepl("se\\(rmean\\)", names(mm))]
  c(est = unname(mm[est_nm][1]), se = unname(mm[se_nm][1]))
}

# Scenario grid: distribution, ties, censoring level, sample size per arm.
scenarios <- expand.grid(
  n    = c(50, 200, 1000),
  dist = c("exp", "weibull"),
  ties = c(FALSE, TRUE),
  cens = c("light", "heavy"),
  stringsAsFactors = FALSE
)

results     <- vector("list", 0L)
results_two <- vector("list", 0L)

for (s in seq_len(nrow(scenarios))) {
  cfg <- scenarios[s, ]
  set.seed(1000 + s)

  n     <- cfg$n
  arm   <- rep(0:1, each = n)
  rate0 <- 0.10
  rate1 <- 0.07

  if (cfg$dist == "exp") {
    surv_t <- c(rexp(n, rate0), rexp(n, rate1))
  } else {
    surv_t <- c(rweibull(n, shape = 1.3, scale = 1 / rate0),
                rweibull(n, shape = 1.3, scale = 1 / rate1))
  }

  cens_rate <- if (cfg$cens == "heavy") 0.10 else 0.03
  cens_t    <- rexp(2 * n, rate = cens_rate)
  obs       <- pmin(surv_t, cens_t)
  status    <- as.integer(surv_t <= cens_t)

  if (cfg$ties) {
    obs <- round(obs, 1)
    obs[obs <= 0] <- 0.1
  }

  tau_max <- min(max(obs[arm == 0]), max(obs[arm == 1]))
  tau     <- 0.8 * tau_max

  # survRM2 reference for both arms (computed once per scenario)
  ref_rm2 <- list(`0` = c(NA, NA, NA, NA), `1` = c(NA, NA, NA, NA))
  ur <- NULL
  if (has_survRM2) {
    invisible(utils::capture.output(
      obj <- survRM2::rmst2(obs, status, arm, tau = tau)
    ))
    ref_rm2[["1"]] <- as.numeric(obj$RMST.arm1$rmst)
    ref_rm2[["0"]] <- as.numeric(obj$RMST.arm0$rmst)
    ur <- obj$unadjusted.result
  }

  for (a in 0:1) {
    idx <- arm == a
    tt  <- obs[idx]
    ss  <- status[idx]

    rf  <- rmst_fast(tt, ss, tau = tau, presorted = FALSE)
    ord <- order(tt)
    rf2 <- rmst_fast(tt[ord], ss[ord], tau = tau, presorted = TRUE)
    rs  <- rmst_survival(tt, ss, tau)
    r2  <- ref_rm2[[as.character(a)]]

    results[[length(results) + 1L]] <- data.frame(
      scenario     = s,
      dist         = cfg$dist,
      ties         = cfg$ties,
      cens         = cfg$cens,
      n_arm        = n,
      arm          = a,
      tau          = tau,
      est_fast     = unname(rf["rmst"]),
      se_fast      = unname(rf["std.err"]),
      presort_diff = max(abs(as.numeric(rf) - as.numeric(rf2))),
      d_est_rm2    = abs(unname(rf["rmst"])    - r2[1]),
      d_se_rm2     = abs(unname(rf["std.err"]) - r2[2]),
      d_lo_rm2     = abs(unname(rf["lower"])   - r2[3]),
      d_hi_rm2     = abs(unname(rf["upper"])   - r2[4]),
      d_est_surv   = abs(unname(rf["rmst"])    - rs["est"]),
      d_se_surv    = abs(unname(rf["std.err"]) - rs["se"]),
      stringsAsFactors = FALSE
    )
  }

  # Two-group mode: contrasts vs survRM2 unadjusted.result
  # ur rows: 1 = "RMST (arm=1)-(arm=0)", 2 = "RMST (arm=1)/(arm=0)"
  # ur cols: 1 = Est., 2 = lower .95, 3 = upper .95, 4 = p
  if (has_survRM2 && !is.null(ur)) {
    tg    <- rmst_fast(obs, status, group = arm, control = 0,
                       tau = tau, presorted = FALSE)
    ord_p <- order(obs)
    tg2   <- rmst_fast(obs[ord_p], status[ord_p], group = arm[ord_p],
                       control = 0, tau = tau, presorted = TRUE)

    results_two[[length(results_two) + 1L]] <- data.frame(
      scenario     = s,
      dist         = cfg$dist,
      ties         = cfg$ties,
      cens         = cfg$cens,
      n_arm        = n,
      tau          = tau,
      presort_diff = max(abs(as.numeric(tg) - as.numeric(tg2))),
      d_rctrl      = abs(unname(tg["rmst.ctrl"])   - ref_rm2[["0"]][1]),
      d_rtrt       = abs(unname(tg["rmst.trt"])    - ref_rm2[["1"]][1]),
      d_diff       = abs(unname(tg["diff"])        - ur[1, 1]),
      d_diff_lo    = abs(unname(tg["diff.lower"])  - ur[1, 2]),
      d_diff_hi    = abs(unname(tg["diff.upper"])  - ur[1, 3]),
      d_diff_p     = abs(unname(tg["p.diff"])      - ur[1, 4]),
      d_ratio      = abs(unname(tg["ratio"])       - ur[2, 1]),
      d_ratio_lo   = abs(unname(tg["ratio.lower"]) - ur[2, 2]),
      d_ratio_hi   = abs(unname(tg["ratio.upper"]) - ur[2, 3]),
      d_ratio_p    = abs(unname(tg["p.ratio"])     - ur[2, 4]),
      stringsAsFactors = FALSE
    )
  }
}

res     <- do.call(rbind, results)
res_two <- if (length(results_two) > 0L) do.call(rbind, results_two) else NULL

# Tolerances
tol_est     <- 1e-6
tol_se      <- 1e-6
tol_presort <- 1e-10

max_presort <- max(res$presort_diff, na.rm = TRUE)
max_est_rm2 <- max(res$d_est_rm2,  na.rm = TRUE)
max_se_rm2  <- max(res$d_se_rm2,   na.rm = TRUE)
max_lo_rm2  <- max(res$d_lo_rm2,   na.rm = TRUE)
max_hi_rm2  <- max(res$d_hi_rm2,   na.rm = TRUE)
max_est_srv <- max(res$d_est_surv, na.rm = TRUE)
max_se_srv  <- max(res$d_se_surv,  na.rm = TRUE)

# Two-group contrast diffs vs survRM2
if (!is.null(res_two)) {
  max_two_presort <- max(res_two$presort_diff, na.rm = TRUE)
  max_two_arm     <- max(res_two$d_rctrl, res_two$d_rtrt, na.rm = TRUE)
  max_two_diff    <- max(res_two$d_diff, res_two$d_diff_lo,
                         res_two$d_diff_hi, na.rm = TRUE)
  max_two_diff_p  <- max(res_two$d_diff_p, na.rm = TRUE)
  max_two_ratio   <- max(res_two$d_ratio, res_two$d_ratio_lo,
                         res_two$d_ratio_hi, na.rm = TRUE)
  max_two_ratio_p <- max(res_two$d_ratio_p, na.rm = TRUE)
} else {
  max_two_presort <- NA_real_; max_two_arm     <- NA_real_
  max_two_diff    <- NA_real_; max_two_diff_p  <- NA_real_
  max_two_ratio   <- NA_real_; max_two_ratio_p <- NA_real_
}

checks <- c(
  presorted_invariant = max_presort < tol_presort,
  est_vs_survRM2      = !has_survRM2 || max_est_rm2 < tol_est,
  se_vs_survRM2       = !has_survRM2 || max_se_rm2  < tol_se,
  ci_vs_survRM2       = !has_survRM2 || (max_lo_rm2 < tol_est && max_hi_rm2 < tol_est),
  est_vs_survival     = is.na(max_est_srv) || max_est_srv < tol_est,
  se_vs_survival      = is.na(max_se_srv)  || max_se_srv  < tol_se,
  twogroup_presorted  = is.na(max_two_presort) || max_two_presort < tol_presort,
  twogroup_per_arm    = is.na(max_two_arm)     || max_two_arm     < tol_est,
  twogroup_diff       = is.na(max_two_diff)    || max_two_diff    < tol_est,
  twogroup_diff_p     = is.na(max_two_diff_p)  || max_two_diff_p  < tol_se,
  twogroup_ratio      = is.na(max_two_ratio)   || max_two_ratio   < tol_est,
  twogroup_ratio_p    = is.na(max_two_ratio_p) || max_two_ratio_p < tol_se
)
all_pass <- all(checks)

# Console + log summary
summary_lines <- c(
  sprintf("Scenarios x arms compared : %d", nrow(res)),
  sprintf("max |presorted diff|      : %.3e (tol %.0e)", max_presort, tol_presort),
  sprintf("max |est - survRM2|       : %.3e (tol %.0e)", max_est_rm2, tol_est),
  sprintf("max |se  - survRM2|       : %.3e (tol %.0e)", max_se_rm2,  tol_se),
  sprintf("max |lower - survRM2|     : %.3e", max_lo_rm2),
  sprintf("max |upper - survRM2|     : %.3e", max_hi_rm2),
  sprintf("max |est - survival|      : %.3e (tol %.0e)", max_est_srv, tol_est),
  sprintf("max |se  - survival|      : %.3e (tol %.0e)", max_se_srv,  tol_se),
  "",
  sprintf("Two-group scenarios       : %d", if (is.null(res_two)) 0L else nrow(res_two)),
  sprintf("max |two-grp presort diff|: %.3e (tol %.0e)", max_two_presort, tol_presort),
  sprintf("max |per-arm - survRM2|   : %.3e (tol %.0e)", max_two_arm,     tol_est),
  sprintf("max |diff/CI - survRM2|   : %.3e (tol %.0e)", max_two_diff,    tol_est),
  sprintf("max |p.diff - survRM2|    : %.3e (tol %.0e)", max_two_diff_p,  tol_se),
  sprintf("max |ratio/CI - survRM2|  : %.3e (tol %.0e)", max_two_ratio,   tol_est),
  sprintf("max |p.ratio - survRM2|   : %.3e (tol %.0e)", max_two_ratio_p, tol_se),
  "",
  "Per-check results:",
  paste(sprintf("  %-20s %s", names(checks), ifelse(checks, "PASS", "FAIL")),
        collapse = "\n"),
  "",
  sprintf("OVERALL: %s", ifelse(all_pass, "PASS", "FAIL"))
)

cat(paste(summary_lines, collapse = "\n"), "\n")
cat(paste(summary_lines, collapse = "\n"), "\n\n", file = log_path, append = TRUE)

# Worst few rows for inspection
ord_worst <- order(-pmax(res$d_est_rm2, res$d_se_rm2))
worst     <- utils::head(res[ord_worst, ], 5L)
cat("Worst 5 rows by |diff| vs survRM2:\n", file = log_path, append = TRUE)
cat(paste(capture.output(print(worst)), collapse = "\n"),
    "\n\n", file = log_path, append = TRUE)

if (!is.null(res_two)) {
  ord_w2 <- order(-pmax(res_two$d_diff, res_two$d_ratio))
  worst2 <- utils::head(res_two[ord_w2, ], 5L)
  cat("Worst 5 two-group rows by |diff|/|ratio| vs survRM2:\n",
      file = log_path, append = TRUE)
  cat(paste(capture.output(print(worst2)), collapse = "\n"),
      "\n\n", file = log_path, append = TRUE)
}

cat("\n--- sessionInfo() ---\n", file = log_path, append = TRUE)
cat(paste(capture.output(sessionInfo()), collapse = "\n"),
    "\n", file = log_path, append = TRUE)

invisible(all_pass)
