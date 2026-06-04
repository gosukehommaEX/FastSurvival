# Cross-validation of ahr_fast() against the archived AHR package
# (Brueckner 2018), which is the reference implementation of the
# Kalbfleisch-Prentice average hazard ratio and the exact package used by
# Dormuth et al. (2024). This is a reproducibility artifact, not a throwaway.
#
# How to run (from the package root, with the package built):
#   devtools::load_all()      # or library(FastSurvival)
#   source("tools/verify_ahr_fast.R")
#
# The archived AHR package is not on CRAN. Install it once from the read-only
# CRAN mirror (it compiles C++), for example:
#   remotes::install_github("cran/AHR")
# If AHR is not installed the script reports that and stops gracefully.

log_path <- file.path("tools", "benchmarks", "verify_ahr_fast_log.txt")
if (!dir.exists(dirname(log_path))) dir.create(dirname(log_path), recursive = TRUE)
cat("ahr_fast vs AHR cross-validation\n", file = log_path, append = FALSE)
cat(format(Sys.time()), "\n\n", file = log_path, append = TRUE)

if (!requireNamespace("AHR", quietly = TRUE) ||
    !requireNamespace("survival", quietly = TRUE)) {
  cat("Package 'AHR' (and/or 'survival') not installed; skipping comparison.\n",
      file = log_path, append = TRUE)
  message("Install AHR with remotes::install_github(\"cran/AHR\") to run this.")
} else {

  # One scenario: returns the maximum absolute differences between ahr_fast and
  # AHR::ahrKM for theta1, theta2, the AHR, Var(theta1), the theta-scale Z (the
  # comparison-group component) and the log(AHR)-scale Z.
  compare_one <- function(seed, ne, rate0, rate1, cens_rate, tau) {
    set.seed(seed)
    o0 <- rexp(ne, rate0)
    o1 <- rexp(ne, rate1)
    cc <- rexp(2 * ne, cens_rate)
    yy <- pmin(c(o0, o1), cc)
    dd <- as.integer(c(o0, o1) <= cc)
    zz <- rep(c(0, 1), each = ne)

    fast <- ahr_fast(yy, dd, zz, tau = tau)

    dat <- data.frame(Y = yy, D = dd, Z = zz)
    ref <- AHR::ahrKM(tau, survival::Surv(Y, D) ~ Z, dat)

    z_loghr <- if (is.null(ref$Z.loghr)) NA_real_ else
      abs(fast$z.loghr - as.numeric(ref$Z.loghr)[1])

    c(theta1 = abs(fast$theta[[1]] - ref$theta[1]),
      theta2 = abs(fast$theta[[2]] - ref$theta[2]),
      ahr = abs(fast$ahr - ref$hr[[1]]),
      var.theta1 = abs(fast$var.theta1 - ref$cov.theta[1, 1]),
      var.theta2 = abs(fast$var.theta2 - ref$cov.theta[2, 2]),
      z.theta = abs(fast$z - ref$Z.theta[2]),
      z.loghr = z_loghr)
  }

  scenarios <- list(
    list(seed = 1,  ne = 100, rate0 = 0.10, rate1 = 0.18, cens_rate = 0.05, tau = 8),
    list(seed = 2,  ne = 200, rate0 = 0.10, rate1 = 0.10, cens_rate = 0.05, tau = 6),
    list(seed = 3,  ne = 150, rate0 = 0.15, rate1 = 0.08, cens_rate = 0.04, tau = 10),
    list(seed = 4,  ne = 300, rate0 = 0.12, rate1 = 0.20, cens_rate = 0.03, tau = 7)
  )

  for (sc in scenarios) {
    d <- do.call(compare_one, sc)
    cat(sprintf("seed=%d ne=%d rate0=%.2f rate1=%.2f tau=%g\n",
                sc$seed, sc$ne, sc$rate0, sc$rate1, sc$tau),
        file = log_path, append = TRUE)
    cat("  max abs diff:\n", file = log_path, append = TRUE)
    for (nm in names(d)) {
      cat(sprintf("    %-12s %.3e\n", nm, d[[nm]]), file = log_path, append = TRUE)
    }
    cat("\n", file = log_path, append = TRUE)
  }

  cat("All differences should be ~1e-8 or smaller (same estimator, ",
      "independent code paths).\n\n", file = log_path, append = TRUE)
}

cat("\n--- sessionInfo() ---\n", file = log_path, append = TRUE)
cat(paste(capture.output(sessionInfo()), collapse = "\n"),
    file = log_path, append = TRUE)
cat("\n", file = log_path, append = TRUE)

message("Wrote ", log_path)
