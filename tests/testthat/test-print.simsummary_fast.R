# ------------------------------------------------------------------ #
#  print.simsummary_fast
# ------------------------------------------------------------------ #

# Two looks, four simulated trials, Z-mode efficacy boundaries.
make_summary <- function() {
  df <- data.frame(sim = rep(1:4, each = 2), look = rep(1:2, 4),
                   logrank.z = c(-2.5, -3.1, -1.0, -1.5, -3.0, -2.2, 0.5, -0.4),
                   n.event = rep(c(50, 100), 4), cutoff = rep(c(12, 24), 4))
  simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
}

test_that("print.simsummary_fast returns its input invisibly", {
  df <- data.frame(sim = 1:4, look = 1L, logrank.z = c(-2.5, -1.0, -3.0, 0.5),
                   n.event = c(50, 52, 48, 55), cutoff = rep(24, 4))
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = -1.96)
  expect_output(print(res), "Group-Sequential Operating Characteristics")
  # capture.output keeps the printed report out of the test log
  utils::capture.output(vis <- withVisible(print(res)))
  expect_false(vis$visible)
  expect_identical(vis$value, res)
})

test_that("print.simsummary_fast prints a selection of rows as a report", {
  res <- make_summary()
  sub <- res[res$look != "overall", ]
  expect_s3_class(sub, "simsummary_fast")
  expect_false(is.null(attr(sub, "boundary")))
  expect_output(print(sub), "Stopping Boundaries: Look by Look")
})

test_that("print.simsummary_fast prints a selection of columns as a data frame", {
  res <- make_summary()
  # Columns of the report are missing.
  sub1 <- res[, c("look", "cum.reject")]
  expect_s3_class(sub1, "simsummary_fast")
  out1 <- utils::capture.output(v <- withVisible(print(sub1)))
  expect_false(any(grepl("Group-Sequential", out1, fixed = TRUE)))
  expect_true(any(grepl("cum.reject", out1, fixed = TRUE)))
  expect_false(v$visible)
  expect_identical(v$value, sub1)
  # All the columns of the report are kept, but not the boundary settings.
  sub2 <- res[, c("population", "look", "cum.reject", "prob.stop.efficacy")]
  expect_null(attr(sub2, "boundary"))
  out2 <- utils::capture.output(print(sub2))
  expect_false(any(grepl("Group-Sequential", out2, fixed = TRUE)))
  expect_true(any(grepl("prob.stop.efficacy", out2, fixed = TRUE)))
})

test_that("print.simsummary_fast prints a selection without look rows as a data frame", {
  res <- make_summary()
  sub <- res[res$look == "overall", ]
  out <- utils::capture.output(print(sub))
  expect_false(any(grepl("Group-Sequential", out, fixed = TRUE)))
  expect_true(any(grepl("overall", out, fixed = TRUE)))
})

test_that("print.simsummary_fast prints a selection of some looks as a data frame", {
  res <- make_summary()
  # Look 2 only: the boundaries refer to both looks of the design.
  out1 <- utils::capture.output(print(res[res$look %in% c("2", "overall"), ]))
  expect_false(any(grepl("Stopping Boundaries", out1, fixed = TRUE)))
  expect_true(any(grepl("cum.reject", out1, fixed = TRUE)))
  # Reordered rows.
  out2 <- utils::capture.output(print(res[c(2, 1, 3), ]))
  expect_false(any(grepl("Stopping Boundaries", out2, fixed = TRUE)))
  # Two of three looks: printed without an error.
  df3 <- data.frame(sim = rep(1:4, each = 3), look = rep(1:3, 4),
                    logrank.z = rep(c(-1, -2, -3), 4),
                    n.event = rep(c(30, 60, 90), 4),
                    cutoff = rep(c(6, 12, 18), 4))
  res3 <- simsummary_fast(df3, eff.col = "logrank.z",
                          efficacy = c(-3.5, -2.8, -1.96))
  out3 <- utils::capture.output(
    v <- withVisible(print(res3[res3$look %in% c("1", "3", "overall"), ])))
  expect_false(v$visible)
  expect_false(any(grepl("Stopping Boundaries", out3, fixed = TRUE)))
})

test_that("print.simsummary_fast pairs each look with its own boundary", {
  res <- make_summary()
  out <- utils::capture.output(print(res))
  # Header of the boundary columns (the statistic need not be a Z-score).
  expect_true(any(grepl("Efficacy Bound", out, fixed = TRUE)))
  tab_line <- function(lk) {
    out[grepl(paste0("^ +", lk, " +[0-9.]+ +[0-9.]+ +-"), out)][1L]
  }
  expect_match(tab_line(1), "-2.8000", fixed = TRUE)
  expect_match(tab_line(2), "-1.9600", fixed = TRUE)
  # A selection of whole blocks (here the single population) keeps the report.
  sub <- res[res$population == res$population[1L], ]
  expect_output(print(sub), "Stopping Boundaries: Look by Look")
})
