# The expected values were computed independently in Python. For the data of
# print_lr_data() the log-rank test (censored times at risk at their own time)
# gives O1 = 5, E1 = 7.4935065, V1 = 1.0066317, a chi-square of 6.1766131 with
# p = 0.0129450, and Z = -2.4852793 with the one-sided p = 0.0064725.
print_lr_data <- function() {
  list(time  = c(1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12),
       event = c(1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1),
       group = c(0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1))
}

test_that("print.survdiff_fast prints the two-sided log-rank test and returns x invisibly", {
  d <- print_lr_data()
  x <- survdiff_fast(d$time, d$event, d$group, control = 0)
  expect_equal(attr(x, "E1"), 7.4935064935, tolerance = 1e-9)
  expect_equal(attr(x, "V1"), 1.0066317160, tolerance = 1e-9)
  expect_output(v <- withVisible(print(x)), "Log-rank test (two-group)",
                fixed = TRUE)
  expect_false(v$visible)
  expect_identical(v$value, x)
  out <- utils::capture.output(print(x))
  expect_true(any(grepl("N = 12,  control = 0", out, fixed = TRUE)))
  expect_true(any(grepl("^treatment +5 +7\\.4935 ", out)))
  expect_true(any(grepl(
    paste0("Chi-square = ", format(signif(6.1766130890, 4)), " on 1 df,  ",
           "p-value = ", format.pval(0.0129449835, digits = 4)),
    out, fixed = TRUE)))
})

test_that("print.survdiff_fast prints the one-sided Z and its lower-tail p-value", {
  d <- print_lr_data()
  x <- survdiff_fast(d$time, d$event, d$group, control = 0, side = 1)
  expect_equal(as.numeric(x), -2.4852792779, tolerance = 1e-9)
  expect_output(print(x),
                paste0("Z = ", format(signif(-2.4852792779, 4)),
                       ",  one-sided p-value = ",
                       format.pval(0.0064724917, digits = 4)),
                fixed = TRUE)
})

test_that("print.survdiff_fast names the weight and shows only the observed counts", {
  d <- print_lr_data()
  x <- survdiff_fast(d$time, d$event, d$group, control = 0, weight = "fh",
                     gamma = 1)
  out <- utils::capture.output(print(x))
  expect_true(any(out == "Fleming-Harrington weighted log-rank test (two-group)"))
  expect_true(any(grepl("Observed", out, fixed = TRUE)))
  expect_false(any(grepl("Expected", out, fixed = TRUE)))
  xs <- survdiff_fast(d$time, d$event, d$group, control = 0, weight = "gehan",
                      strata = rep(1:2, each = 6))
  expect_output(print(xs), paste0("Stratified Gehan-Breslow weighted log-rank ",
                                  "test (two-group, 2 strata)"), fixed = TRUE)
})

test_that("print labels the unweighted stratified log-rank test", {
  skip_if_not_installed("survival")
  ov  <- survival::ovarian
  fit <- survdiff_fast(ov$futime, ov$fustat, ov$rx, 1, strata = ov$resid.ds)
  expect_output(print(fit), "Stratified log-rank test")
})

test_that("print.survdiff_fast reports an unavailable statistic", {
  d <- print_lr_data()
  x <- survdiff_fast(d$time, d$event, d$group, control = 0)
  x_na <- x
  x_na[] <- NA_real_
  expect_output(print(x_na), "Test statistic not available (V1 = 0).",
                fixed = TRUE)
})
