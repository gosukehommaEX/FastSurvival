make_print_data <- function() {
  simdata_fast(nsim = 1, n = c(60, 60), a.time = c(0, 6), a.prop = 1,
               e.median = list(8, 14), d.hazard = 0.02, seed = 501)
}

test_that("print.coxph_fast labels the p-value by the side of the test", {
  d <- make_print_data()
  x2 <- coxph_fast(d$tte, d$event, d$group, control = 1)
  expect_output(v <- withVisible(print(x2)), "Pr(>|z|)", fixed = TRUE)
  expect_output(print(x2), "Cox proportional hazards model for two groups",
                fixed = TRUE)
  expect_output(print(x2), "ties = efron", fixed = TRUE)
  expect_false(v$visible)
  expect_identical(v$value, x2)
  x1 <- coxph_fast(d$tte, d$event, d$group, control = 1, side = 1)
  out1 <- utils::capture.output(print(x1))
  expect_true(any(grepl("Pr(<z)", out1, fixed = TRUE)))
  expect_false(any(grepl("Pr(>|z|)", out1, fixed = TRUE)))
  expect_true(any(grepl("one.sided", out1, fixed = TRUE)))
})

test_that("print.coxph_fast explains an unavailable estimate", {
  # All events are in the control group, so the estimate is infinite
  na_x <- coxph_fast(1:6, c(1, 1, 1, 0, 0, 0), c(0, 0, 0, 1, 1, 1),
                     control = 0)
  expect_output(print(na_x), "no events, or an infinite estimate",
                fixed = TRUE)
  nc_x <- structure(na_x, converged = FALSE)
  expect_output(print(nc_x), "the iteration did not converge", fixed = TRUE)
})
