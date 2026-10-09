# Scaling of the FastSurvival pipeline with the number of simulated trials,
# and reproducibility of batched execution with dqrng streams.
#
# Part 1 times simdata_fast() and analysis_fast() (two event-driven looks; the
# log-rank and RMST statistics and the max-combo test timed separately) for an
# increasing number of simulated trials and records the size of the simulated
# data. The max-combo test is timed twice: with every p-value integrated, and
# with mc.alpha set to the nominal levels of a group-sequential design, which
# integrates only the p-values that the Bonferroni bounds do not decide. Part 2
# runs the same study in batches, once sequentially and once on a parallel
# cluster, and checks that the combined results are identical.
# Results are written to tools/paper/output/.
#
# Run from the package root after installing FastSurvival 1.1.0 from CRAN
# (checked by machine_info.R):
#   source("tools/paper/data/bench_scaling.R")

library(FastSurvival)
library(parallel)
source(file.path("tools", "paper", "data", "machine_info.R"))

out_dir <- file.path("tools", "paper", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

design <- list(n = c(300, 300), a.time = c(0, 12), a.rate = 600 / 12,
               e.hazard = list(c(0.06, 0.06), c(0.06, 0.035)),
               e.time = c(0, 4, Inf), d.hazard = 0.001)
looks <- c(200, 350)

# Nominal one-sided levels of a Lan-DeMets O'Brien-Fleming design at these
# looks, used for the max-combo shortcut (mc.alpha) in Part 1.
mc_alpha <- stats::pnorm(-gsDesign::gsDesign(
  k = 2, test.type = 1, alpha = 0.025, timing = looks / max(looks),
  sfu = gsDesign::sfLDOF)$upper$bound)

run_once <- function(nsim, seed, stream = NULL) {
  d <- do.call(simdata_fast, c(list(nsim = nsim, seed = seed, stream = stream),
                               design))
  analysis_fast(d, control = 1, event.looks = looks,
                stat = c("logrank", "maxcombo", "rmst"), tau = 18, side = 1)
}

# ---- Part 1: time and memory against nsim ----------------------------------
# The deterministic statistics (log-rank and RMST, computed in the fused C++
# loop) and the max-combo test (whose four-weight p-value is a multivariate
# normal integral evaluated by mvtnorm in R, one call per simulated trial and
# look) are timed separately.
grid <- c(1000, 2500, 5000, 10000)
part1 <- do.call(rbind, lapply(grid, function(ns) {
  t_gen <- system.time(
    d <- do.call(simdata_fast, c(list(nsim = ns, seed = 1), design))
  )[["elapsed"]]
  t_det <- system.time(
    analysis_fast(d, control = 1, event.looks = looks,
                  stat = c("logrank", "rmst"), tau = 18, side = 1)
  )[["elapsed"]]
  set.seed(1)
  t_mc <- system.time(
    r_mc <- analysis_fast(d, control = 1, event.looks = looks,
                          stat = "maxcombo", side = 1)
  )[["elapsed"]]
  set.seed(1)
  t_mb <- system.time(
    r_mb <- analysis_fast(d, control = 1, event.looks = looks,
                          stat = "maxcombo", side = 1, mc.alpha = mc_alpha)
  )[["elapsed"]]
  a_row <- mc_alpha[r_mc$look]
  data.frame(nsim = ns, rows = nrow(d),
             data_mb = as.numeric(object.size(d)) / 2^20,
             generate_sec = t_gen, logrank_rmst_sec = t_det,
             maxcombo_sec = t_mc, maxcombo_mc_alpha_sec = t_mb,
             share_integrated = mean(r_mb$maxcombo.p.exact, na.rm = TRUE),
             decisions_differ = sum((r_mb$maxcombo.p <= a_row) !=
                                      (r_mc$maxcombo.p <= a_row),
                                    na.rm = TRUE))
}))
print(part1)
write.csv(part1, file.path(out_dir, "bench_scaling.csv"), row.names = FALSE)

# ---- Part 2: batches from independent streams -------------------------------
# Each batch uses dqrng stream b for the data and set.seed(b) for the
# GenzBretz integration of the max-combo p-value, so that every column is
# reproducible whatever worker runs the batch.
n_batch   <- 8
per_batch <- 2500
run_batch <- function(b) {
  set.seed(b)
  r <- run_once(per_batch, seed = 2026, stream = b)
  r$sim <- r$sim + (b - 1) * per_batch
  r
}
t_seq <- system.time(
  seq_res <- do.call(rbind, lapply(seq_len(n_batch), run_batch))
)[["elapsed"]]

n_cores <- max(1L, min(4L, detectCores() - 1L))
cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(FastSurvival))
clusterExport(cl, c("run_once", "design", "looks", "per_batch"))
t_par <- system.time(
  par_res <- do.call(rbind, parLapply(cl, rev(seq_len(n_batch)), run_batch))
)[["elapsed"]]
stopCluster(cl)
par_res <- par_res[order(par_res$sim, par_res$look), ]
rownames(par_res) <- NULL
rownames(seq_res) <- NULL

num_cols <- names(seq_res)[vapply(seq_res, is.numeric, logical(1))]
col_diff <- vapply(num_cols, function(cn) {
  max(abs(seq_res[[cn]] - par_res[[cn]]), na.rm = TRUE)
}, numeric(1))
part2 <- data.frame(
  batches = n_batch, trials = n_batch * per_batch, cores = n_cores,
  sequential_sec = t_seq, parallel_sec = t_par,
  identical_results = identical(seq_res, par_res),
  max_abs_difference = max(col_diff)
)
print(part2)
write.csv(part2, file.path(out_dir, "bench_batches.csv"), row.names = FALSE)
write.csv(data.frame(column = num_cols, max_abs_difference = col_diff),
          file.path(out_dir, "bench_batches_columns.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "bench_scaling_sessionInfo.txt"))
