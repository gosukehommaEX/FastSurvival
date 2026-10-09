# Check of the mc.alpha shortcut of analysis_fast() in bench_scaling.R.
#
# In bench_scaling.R the decisions at the nominal levels differed in a few
# analyses between the runs with and without mc.alpha. The four-weight
# max-combo p-value is integrated by the randomized GenzBretz algorithm, and
# the two runs integrate a given analysis with different random numbers. For
# 10,000 simulated trials of the bench_scaling.R design, this script compares
# the decisions of (1) all p-values integrated with set.seed(1), (2) the run
# with mc.alpha and set.seed(1), and (3) all p-values integrated again with
# set.seed(2), and integrates the p-values of the analyses whose decisions
# differ between (1) and (2) once more with a smaller error tolerance.
# Results are written to tools/paper/output/check_mc_alpha_summary.csv and
# check_mc_alpha_rows.csv.
#
# Run from the package root after installing FastSurvival 1.1.0 from CRAN
# (checked by machine_info.R):
#   source("tools/paper/data/check_mc_alpha.R")

library(FastSurvival)
source(file.path("tools", "paper", "data", "machine_info.R"))

out_dir <- file.path("tools", "paper", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Design, looks, and nominal levels as in bench_scaling.R.
design <- list(n = c(300, 300), a.time = c(0, 12), a.rate = 600 / 12,
               e.hazard = list(c(0.06, 0.06), c(0.06, 0.035)),
               e.time = c(0, 4, Inf), d.hazard = 0.001)
looks <- c(200, 350)
mc_alpha <- stats::pnorm(-gsDesign::gsDesign(
  k = 2, test.type = 1, alpha = 0.025, timing = looks / max(looks),
  sfu = gsDesign::sfLDOF)$upper$bound)

nsim <- 10000
d <- do.call(simdata_fast, c(list(nsim = nsim, seed = 1), design))
set.seed(1)
r1 <- analysis_fast(d, control = 1, event.looks = looks, stat = "maxcombo",
                    side = 1)
set.seed(1)
rb <- analysis_fast(d, control = 1, event.looks = looks, stat = "maxcombo",
                    side = 1, mc.alpha = mc_alpha)
set.seed(2)
r2 <- analysis_fast(d, control = 1, event.looks = looks, stat = "maxcombo",
                    side = 1)
a <- mc_alpha[r1$look]
dec1 <- r1$maxcombo.p <= a
decb <- rb$maxcombo.p <= a
dec2 <- r2$maxcombo.p <= a
differ <- which(dec1 != decb)

# Integrate the analyses whose decisions differ again with a smaller
# tolerance and more points.
sims <- unique(r1$sim[differ])
set.seed(3)
rp <- analysis_fast(d[d$sim %in% sims, ], control = 1, event.looks = looks,
                    stat = "maxcombo", side = 1, abseps = 1e-7, maxpts = 1e6)
key  <- paste(r1$sim, r1$look)
keyp <- paste(rp$sim, rp$look)
rows <- data.frame(
  sim = r1$sim[differ], look = r1$look[differ], level = a[differ],
  p_seed1 = r1$maxcombo.p[differ], p_seed2 = r2$maxcombo.p[differ],
  p_mc_alpha = rb$maxcombo.p[differ],
  integrated = rb$maxcombo.p.exact[differ],
  p_precise = rp$maxcombo.p[match(key[differ], keyp)]
)
noise <- abs(r1$maxcombo.p - r2$maxcombo.p)
summ <- data.frame(
  nsim = nsim, analyses = nrow(r1),
  integrated_mc_alpha = sum(rb$maxcombo.p.exact, na.rm = TRUE),
  decisions_differ_mc_alpha = length(differ),
  decisions_differ_seed = sum(dec1 != dec2, na.rm = TRUE),
  max_gap_differ = if (length(differ)) {
    max(abs(rows$p_seed1 - rows$level))
  } else {
    0
  },
  max_abs_noise = max(noise, na.rm = TRUE),
  q99_abs_noise = unname(stats::quantile(noise, 0.99, na.rm = TRUE)),
  precise_agrees_mc_alpha = sum((rows$p_precise <= rows$level) ==
                                  (rows$p_mc_alpha <= rows$level)),
  precise_agrees_seed1 = sum((rows$p_precise <= rows$level) ==
                               (rows$p_seed1 <= rows$level))
)
print(summ)
print(rows)
write.csv(summ, file.path(out_dir, "check_mc_alpha_summary.csv"),
          row.names = FALSE)
write.csv(rows, file.path(out_dir, "check_mc_alpha_rows.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "check_mc_alpha_sessionInfo.txt"))
