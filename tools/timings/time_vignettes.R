# Measure the rendering time of every package vignette before CRAN submission.
# Each vignette is rendered in an isolated R process so that timings are not
# contaminated by a warm package load, knitr global state, or one vignette's
# objects leaking into the next. The elapsed time is what matters for CRAN,
# which renders from scratch with no cache.
#
# This is a one-shot measurement, not part of the checked package surface. Run
# it from the package root (the directory that contains the "vignettes" folder)
# with FastSurvival installed.
#
# Set 'fresh' to TRUE to clear each vignette's knitr cache before timing, which
# reproduces the cold render that CRAN performs. Leave it FALSE to time the
# vignettes as they currently stand (warm cache, faster, not CRAN-representative).

fresh <- FALSE

vignette_dir <- "vignettes"
out_dir      <- file.path("tools", "timings")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
log_path <- file.path(out_dir, "time_vignettes_log.txt")

cat("Vignette rendering times\n",
    "Run at: ", format(Sys.time()), "\n",
    "fresh (cold render) = ", fresh, "\n\n",
    file = log_path, append = FALSE, sep = "")

if (!dir.exists(vignette_dir)) {
  stop("Directory '", vignette_dir, "' not found. Run from the package root.")
}

rmd_files <- sort(list.files(vignette_dir, pattern = "\\.Rmd$",
                             full.names = TRUE))
if (length(rmd_files) == 0L) {
  stop("No .Rmd vignettes found in '", vignette_dir, "'.")
}

cat("Found", length(rmd_files), "vignette(s):\n")
cat(paste0("  ", basename(rmd_files), collapse = "\n"), "\n\n")

has_callr <- requireNamespace("callr", quietly = TRUE)
if (!has_callr) {
  warning("Package 'callr' is not installed; rendering in-process instead. ",
          "Timings will be less isolated. Install 'callr' for clean numbers.")
}

# Remove the knitr cache and generated files for one vignette so the next render
# is cold. Only the conventional sidecar locations are touched.
clear_artifacts <- function(rmd) {
  stem <- tools::file_path_sans_ext(rmd)
  targets <- c(file.path(vignette_dir, "cache"),
               paste0(stem, "_cache"),
               paste0(stem, "_files"))
  for (t in targets) {
    if (dir.exists(t)) unlink(t, recursive = TRUE, force = TRUE)
  }
}

# Render one vignette and return its elapsed seconds. The render happens in a
# throwaway output file and a fresh environment; resources resolve relative to
# the vignette's own directory, as during package build.
render_one <- function(rmd) {
  if (has_callr) {
    callr::r(
      function(input) {
        out <- tempfile(fileext = ".html")
        on.exit(unlink(out), add = TRUE)
        tt <- system.time(
          rmarkdown::render(input, output_file = out,
                            quiet = TRUE, envir = new.env())
        )
        as.numeric(tt[["elapsed"]])
      },
      args = list(input = normalizePath(rmd)),
      wd = getwd()
    )
  } else {
    out <- tempfile(fileext = ".html")
    on.exit(unlink(out), add = TRUE)
    tt <- system.time(
      rmarkdown::render(normalizePath(rmd), output_file = out,
                        quiet = TRUE, envir = new.env())
    )
    as.numeric(tt[["elapsed"]])
  }
}

elapsed <- rep(NA_real_, length(rmd_files))
status  <- rep("ok", length(rmd_files))

for (i in seq_along(rmd_files)) {
  rmd <- rmd_files[i]
  if (fresh) clear_artifacts(rmd)
  cat(sprintf("Rendering %-34s ... ", basename(rmd)))
  res <- tryCatch(render_one(rmd), error = function(e) e)
  if (inherits(res, "error")) {
    status[i]  <- "error"
    cat("ERROR\n")
    cat(sprintf("[ERROR] %s : %s\n", basename(rmd), conditionMessage(res)),
        file = log_path, append = TRUE)
  } else {
    elapsed[i] <- res
    cat(sprintf("%7.1f s\n", res))
    cat(sprintf("%-34s %8.1f s\n", basename(rmd), res),
        file = log_path, append = TRUE)
  }
}

# A soft flag for vignettes that may be too slow for comfortable CRAN checks.
# CRAN dislikes long-running vignettes; the usual remedy is precomputation.
slow_threshold <- 30

results <- data.frame(
  vignette = basename(rmd_files),
  status   = status,
  elapsed_s = round(elapsed, 1),
  slow      = !is.na(elapsed) & elapsed > slow_threshold,
  stringsAsFactors = FALSE
)
results <- results[order(-results$elapsed_s, na.last = TRUE), ]
rownames(results) <- NULL

total_s <- sum(elapsed, na.rm = TRUE)
cat("\n", file = log_path, append = TRUE)
cat(sprintf("Total elapsed over %d vignette(s): %.1f s\n",
            sum(status == "ok"), total_s),
    file = log_path, append = TRUE)

cat("\n")
print(results)
cat(sprintf("\nTotal elapsed over %d vignette(s): %.1f s\n",
            sum(status == "ok"), total_s))
if (any(results$slow, na.rm = TRUE)) {
  cat("Flagged as slow (> ", slow_threshold, " s): ",
      paste(results$vignette[which(results$slow)], collapse = ", "), "\n",
      sep = "")
}

cat("\nsessionInfo():\n", file = log_path, append = TRUE)
capture.output(sessionInfo(), file = log_path, append = TRUE)
