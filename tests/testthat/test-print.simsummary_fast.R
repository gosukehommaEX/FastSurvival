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
