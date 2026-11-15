# The expected values were computed independently in Python. With uncensored
# times 1, 3, ..., 19 (control) and 2, 4, ..., 20 (treatment) the
# Kaplan-Meier estimate is 0.5 between the fifth and the sixth event time, so
# the medians are the midpoints 10 and 11 (the convention of
# survival::survfit), and the Cox partial likelihood estimate (no ties) is
# -0.2265158, a hazard ratio of 0.797.

test_that("print.kmcurve_fast prints the two-group summary and returns x invisibly", {
  x <- kmcurve_fast(c(seq(1, 19, by = 2), seq(2, 20, by = 2)), rep(1L, 20),
                    rep(c(0, 1), each = 10), control = 0)
  expect_output(v <- withVisible(print(x)),
                "A kmcurve_fast object (two groups)", fixed = TRUE)
  expect_false(v$visible)
  expect_identical(v$value, x)
  out <- utils::capture.output(print(x))
  expect_true(any(grepl("^ *0 +control +10 +10 +10$", out)))
  expect_true(any(grepl("^ *1 +treatment +10 +10 +11$", out)))
  expect_true(any(out == "Cox PH hazard ratio (treatment vs control): 0.797"))
})

test_that("print.kmcurve_fast prints the one-group summary without a hazard ratio", {
  # Uncensored times 1, ..., 10: the median is the midpoint 5.5.
  x <- kmcurve_fast(1:10, rep(1L, 10))
  out <- utils::capture.output(v <- withVisible(print(x)))
  expect_false(v$visible)
  expect_identical(v$value, x)
  expect_true(any(out == "A kmcurve_fast object (one group)"))
  expect_true(any(grepl("^ *All +single +10 +10 +5.5$", out)))
  expect_false(any(grepl("hazard ratio", out, fixed = TRUE)))
})

test_that("print.kmcurve_fast omits the hazard ratio when the estimate is infinite", {
  # All events are in the control group (times 1 to 5), so the Cox estimate
  # is infinite, the treatment median is not reached, and the control median
  # is the event time 3 (the curve falls from 0.6 to 0.4 there).
  x <- kmcurve_fast(1:10, rep(c(1, 0), each = 5), rep(c(0, 1), each = 5),
                    control = 0)
  out <- utils::capture.output(print(x))
  expect_true(any(grepl("^ *0 +control +5 +5 +3$", out)))
  expect_true(any(grepl("^ *1 +treatment +5 +0 +NA$", out)))
  expect_false(any(grepl("hazard ratio", out, fixed = TRUE)))
})
