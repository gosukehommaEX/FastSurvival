# Write the table of speed gains in vignettes/speed-comparison.Rmd, and the
# date, R version, and operating system in its sentence "measured on ...",
# from the output of tools/benchmark_speed.R (tools/benchmarks/
# speed_results.csv, speed_results_session.txt, and speed_results.log). Only
# the third column of the table and the date, R version, and operating
# system of that sentence are changed; the wording of the first two columns
# is kept.
#
# The speed gain is the column speedup of speed_results.csv (the median time
# of the reference over that of FastSurvival). It is rounded to a multiple of
# 5 below 100, of 50 below 1000, and of 500 from 1000 on, with halves rounded
# up. The rule is first checked against the table of the run of 2026-09-29,
# which it reproduces and which rounding to a multiple of 10 does not.
#
# Run from the package root after tools/benchmark_speed.R:
#   source("tools/update_speed_table.R")
# A summary is written to tools/benchmarks/update_speed_table.log.

gain_label <- function(x) {
  step <- if (x < 100) 5 else if (x < 1000) 50 else 500
  sprintf("~%dx", as.integer(floor(x / step + 0.5) * step))
}

# The run of 2026-09-29 (speed_results.csv) and the table made from it
old_run <- c(38, 26.3, 35.7, 36.6, 298.7, 875, 18.3, 37.2, 357.6, 77.4, 450.2,
             198.1)
old_table <- c("~40x", "~25x", "~35x", "~35x", "~300x", "~900x", "~20x",
               "~35x", "~350x", "~75x", "~450x", "~200x")
by_ten <- sprintf("~%dx", as.integer(floor(old_run / 10 + 0.5) * 10))
stopifnot(identical(vapply(old_run, gain_label, ""), old_table),
          !identical(by_ten, old_table))

bench <- file.path("tools", "benchmarks")
res <- utils::read.csv(file.path(bench, "speed_results.csv"),
                       stringsAsFactors = FALSE)
ses <- readLines(file.path(bench, "speed_results_session.txt"),
                 encoding = "UTF-8")
session <- stats::setNames(sub("^[^:]+: ", "", ses), sub(": .*$", "", ses))
log_txt <- readLines(file.path(bench, "speed_results.log"), encoding = "UTF-8")
stopifnot(any(grepl("n = 500, events = 500,", log_txt, fixed = TRUE)),
          any(grepl("replicates = 1000", log_txt, fixed = TRUE)))
skip_pat <- "^Skipped (\\S+) benchmark \\((\\S+) not installed\\)"
skip <- regmatches(log_txt, regexec(skip_pat, log_txt, perl = TRUE))
skip <- vapply(skip[lengths(skip) > 0L],
               function(m) paste(m[2L], m[3L], sep = "/"), "")

r_ver <- sub("^R version ([0-9]+\\.[0-9]+\\.[0-9]+) .*$", "\\1", session[["R"]])
os_label <- session[["running"]]
if (grepl("^Windows [0-9]+ x64", os_label)) {
  os_label <- sub("^Windows ([0-9]+) x64.*$", "Windows \\1 (x86_64)", os_label)
}

vig <- file.path("vignettes", "speed-comparison.Rmd")
raw <- readBin(vig, "raw", file.size(vig))
stopifnot(!any(raw == as.raw(13L)))
md5_before <- unname(tools::md5sum(vig))
txt <- readLines(vig, encoding = "UTF-8", warn = FALSE)

pat <- paste0("measured on [0-9]{4}-[0-9]{2}-[0-9]{2} with R ",
              "[0-9]+\\.[0-9]+\\.[0-9]+ on .+? for the data set of")
hit <- grep(pat, txt, perl = TRUE)
stopifnot(length(hit) == 1L)
m <- regexpr(pat, txt[hit], perl = TRUE)
regmatches(txt[hit], m) <- paste0("measured on ", session[["date"]], " with R ",
                                  r_ver, " on ", os_label,
                                  " for the data set of")
whole <- paste(txt, collapse = "\n")
stopifnot(grepl(paste("The reported figures are median times of 1,000",
                      "microbenchmark replicates,"), whole, fixed = TRUE),
          grepl("500 subjects below, all of whom have an event.", whole,
                fixed = TRUE))

row_pat <- "^\\| `(.+?)` \\| (.+) \\| (~[0-9]+x) \\|$"
rows <- grep(row_pat, txt, perl = TRUE)
keys <- sub(row_pat, "\\1", txt[rows], perl = TRUE)
stopifnot(setequal(keys, res$method), !anyDuplicated(keys))
changes <- character(0)
for (k in seq_along(rows)) {
  i <- rows[k]
  gain <- res$speedup[res$method == keys[k]]
  old <- sub(row_pat, "\\3", txt[i], perl = TRUE)
  new <- gain_label(gain)
  txt[i] <- paste0("| `", keys[k], "` | ", sub(row_pat, "\\2", txt[i],
                                               perl = TRUE), " | ", new, " |")
  changes <- c(changes, sprintf("  %s: speed gain %s -> %s (was %s)", keys[k],
                                format(gain), new, old))
}

missing <- paste("The `wkm_fast()` row is missing because",
                 "nphsim was not installed when the table was produced.",
                 sep = "\n")
if (identical(skip, "wkm_fast/nphsim")) {
  stopifnot(lengths(regmatches(paste(txt, collapse = "\n"),
                               gregexpr(missing, paste(txt, collapse = "\n"),
                                        fixed = TRUE))) == 1L)
} else {
  stop("the sentence on missing rows must be rewritten for: ",
       paste(skip, collapse = ", "))
}

con <- file(vig, open = "wb")
writeLines(enc2utf8(txt), con, sep = "\n", useBytes = TRUE)
close(con)
md5_after <- unname(tools::md5sum(vig))

summary_lines <- c(
  paste0("update_speed_table.R run on ", format(Sys.time()), " with ",
         R.version.string),
  sprintf(paste("speed_results.csv: %d rows; run of %s, R %s, %s,",
                "FastSurvival %s; skipped %s"),
          nrow(res), session[["date"]], r_ver, os_label,
          session[["FastSurvival"]], paste(skip, collapse = ", ")),
  changes,
  sprintf("vignettes/speed-comparison.Rmd %s -> %s", md5_before, md5_after))
writeLines(summary_lines, file.path(bench, "update_speed_table.log"))
cat(summary_lines, sep = "\n")
