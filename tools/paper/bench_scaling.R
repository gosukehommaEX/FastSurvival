# Scaling of the FastSurvival pipeline with the number of simulated trials,
# and reproducibility of batched execution with dqrng streams.
#
# Part 1 times simdata_fast() + analysis_fast() (two event-driven looks, three
# statistics) for an increasing number of simulated trials and records the
# size of the simulated data. Part 2 runs the same study in batches, once
# sequentially and once on a parallel cluster, and checks that the combined
# results are identical. Results are written to tools/paper/output/.
#
# Run from the package root after installing the package:
#   source("tools/paper/bench_scaling.R")

library(FastSurvival)
library(parallel)

out_dir <- file.path("tools", "paper", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

design <- list(n = c(300, 300), a.time = c(0, 12), a.rate = 600 / 12,
               e.hazard = list(c(0.06, 0.06), c(0.06, 0.035)),
               e.time = c(0, 4, Inf), d.hazard = 0.001)
looks <- c(200, 350)

run_once <- function(nsim, seed, stream = NULL) {
  d <- do.call(simdata_fast, c(list(nsim = nsim, seed = seed, stream = stream),
                               design))
  analysis_fast(d, control = 1, event.looks = looks,
                stat = c("logrank", "maxcombo", "rmst"), tau = 18, side = 1)
}

# ---- Part 1: time and memory against nsim ----------------------------------
grid <- c(1000, 2500, 5000, 10000)
part1 <- do.call(rbind, lapply(grid, function(ns) {
  t_gen <- system.time(
    d <- do.call(simdata_fast, c(list(nsim = ns, seed = 1), design))
  )[["elapsed"]]
  t_ana <- system.time(
    r <- analysis_fast(d, control = 1, event.looks = looks,
                       stat = c("logrank", "maxcombo", "rmst"), tau = 18,
                       side = 1)
  )[["elapsed"]]
  data.frame(nsim = ns, rows = nrow(d),
             data_mb = as.numeric(object.size(d)) / 2^20,
             generate_sec = t_gen, analyze_sec = t_ana,
             total_sec = t_gen + t_ana)
}))
print(part1)
write.csv(part1, file.path(out_dir, "bench_scaling.csv"), row.names = FALSE)

# ---- Part 2: batches from independent streams -------------------------------
n_batch   <- 8
per_batch <- 2500
run_batch <- function(b) {
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

part2 <- data.frame(
  batches = n_batch, trials = n_batch * per_batch, cores = n_cores,
  sequential_sec = t_seq, parallel_sec = t_par,
  identical_results = isTRUE(all.equal(seq_res, par_res))
)
print(part2)
write.csv(part2, file.path(out_dir, "bench_batches.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "bench_scaling_sessionInfo.txt"))
