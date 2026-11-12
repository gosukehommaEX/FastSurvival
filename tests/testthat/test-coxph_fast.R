test_that("coxph_fast returns named vector of length 5", {
  set.seed(1)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_length(res, 5L)
  expect_named(res, c("coef", "exp(coef)", "se(coef)", "lower .95", "upper .95"))
})

test_that("coxph_fast exp(coef) is positive", {
  set.seed(2)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_gt(unname(res["exp(coef)"]), 0)
})

test_that("coxph_fast CI lower <= exp(coef) <= upper", {
  set.seed(3)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_lte(unname(res["lower .95"]), unname(res["exp(coef)"]) + 1e-10)
  expect_gte(unname(res["upper .95"]), unname(res["exp(coef)"]) - 1e-10)
})

test_that("coxph_fast coef = log(exp(coef))", {
  set.seed(4)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_equal(unname(res["coef"]), log(unname(res["exp(coef)"])),
               tolerance = 1e-12)
})

test_that("coxph_fast presorted=TRUE and presorted=FALSE agree", {
  set.seed(5)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)
  ord   <- order(time)

  res_pre <- coxph_fast(time[ord], event[ord], group[ord],
                        control = 1, presorted = TRUE)
  res_uns <- coxph_fast(time, event, group,
                        control = 1, presorted = FALSE)
  expect_equal(unname(res_pre), unname(res_uns), tolerance = 1e-12)
})

test_that("coxph_fast conf.level argument changes CI label names", {
  set.seed(6)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res90 <- coxph_fast(time, event, group, control = 1, conf.level = 0.90)
  expect_named(res90, c("coef", "exp(coef)", "se(coef)", "lower .90", "upper .90"))
})

test_that("coxph_fast handles factor group correctly", {
  set.seed(7)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- factor(rep(c("control", "treatment"), each = n / 2))

  res <- coxph_fast(time, event, group, control = "control")
  expect_gt(unname(res["exp(coef)"]), 0)
  expect_false(any(is.na(res)))
})

test_that("coxph_fast agrees with survival::coxph() (estimate and SE)", {
  skip_if_not_installed("survival")
  ctrl <- survival::coxph.control(timefix = FALSE)
  lung <- survival::lung
  sets <- list(
    ovarian = list(t = survival::ovarian$futime, e = survival::ovarian$fustat,
                   g = survival::ovarian$rx),
    veteran = list(t = survival::veteran$time, e = survival::veteran$status,
                   g = survival::veteran$trt),
    lung = list(t = lung$time, e = as.integer(lung$status == 2), g = lung$sex),
    # Times in months: many tied event times
    lung_month = list(t = ceiling(lung$time / 30.4375),
                      e = as.integer(lung$status == 2), g = lung$sex)
  )
  for (nm in names(sets)) {
    d <- sets[[nm]]
    x <- as.numeric(d$g != min(d$g))
    for (tie in c("efron", "breslow")) {
      fast <- coxph_fast(d$t, d$e, d$g, control = min(d$g), ties = tie)
      ref <- survival::coxph(survival::Surv(d$t, d$e) ~ x, ties = tie,
                             control = ctrl)
      lab <- paste(nm, tie)
      expect_equal(unname(fast["coef"]), unname(stats::coef(ref)),
                   tolerance = 1e-10, label = lab)
      expect_equal(unname(fast["se(coef)"]), sqrt(ref$var[1, 1]),
                   tolerance = 1e-10, label = lab)
      expect_identical(attr(fast, "iter"), as.integer(ref$iter), label = lab)
      expect_true(attr(fast, "converged"), label = lab)
    }
  }
})

test_that("coxph_fast reproduces Therneau and Grambsch (2000), Appendix E", {
  # Test data 1 (p. 323). Published (p. 324, Breslow): beta 1.475285,
  # LL(0) -4.564348, LL(beta) -3.824750, I(beta) 0.634168, iterates 0, 8/5,
  # 1.4727236, 1.4752838, 1.4752849; (p. 327, Efron): LL(0) -4.276666. The
  # expected values to 1e-12 are the exact Newton iterates of the book's own
  # U and I and the fit of Planning/python_checks/cox_exact.py. The printed
  # iterates 2 and 3 and the Efron LL(beta), I(beta), and beta of p. 327
  # differ from the book's own formulas (see cox_exact_output.txt), so they
  # are not tested.
  tt <- c(1, 1, 6, 6, 8, 9)
  ee <- c(1L, 0L, 1L, 1L, 0L, 1L)
  xx <- c(1L, 1L, 1L, 0L, 0L, 0L)
  br <- cox_core(tt, ee, xx, integer(0), efron = FALSE)
  expect_identical(br[c(5, 6)], c(4, 0))
  expect_lt(abs(br[1] - 1.475285), 5e-7)
  expect_lt(abs(br[3] - (-4.564348)), 5e-7)
  expect_lt(abs(br[4] - (-3.824750)), 5e-7)
  expect_lt(abs(1 / br[2]^2 - 0.634168), 5e-7)
  expect_equal(br[1:4], c(1.475284914828868, 1.2557343930105218,
                          -4.564348191467836, -3.8247495050028655),
               tolerance = 1e-12)
  iter_br <- c(1.6, 1.472723532245478, 1.475283960697332, 1.4752849148288683)
  ef <- cox_core(tt, ee, xx, integer(0), efron = TRUE)
  expect_identical(ef[c(5, 6)], c(4, 0))
  expect_lt(abs(ef[3] - (-4.276666)), 5e-7)
  expect_equal(ef[1:4], c(1.6768574855882044, 1.277615576278031,
                          -4.276666119016055, -3.3589748402633477),
               tolerance = 1e-12)
  iter_ef <- c(1.8795180722891567, 1.6701860496115566, 1.6768515659499237,
               1.6768574855882044)
  for (k in 1:4) {
    expect_equal(cox_core(tt, ee, xx, integer(0), FALSE, iter_max = k)[1],
                 iter_br[k], tolerance = 1e-12)
    expect_equal(cox_core(tt, ee, xx, integer(0), TRUE, iter_max = k)[1],
                 iter_ef[k], tolerance = 1e-12)
  }
  expect_identical(cox_core(tt, ee, xx, integer(0), FALSE, iter_max = 1)[1],
                   8 / 5)
  fit <- coxph_fast(tt, ee, xx, control = 0, ties = "breslow")
  expect_identical(unname(fit["coef"]), br[1])
  expect_identical(unname(coxph_fast(tt, ee, xx, control = 0)["coef"]), ef[1])
})

test_that("coxph_fast returns NA when the estimate is infinite", {
  # The finiteness condition L < O_T < H of ?coxph_fast fails in the first
  # three data sets (no root of the score; Planning/python_checks/
  # cox_exact.py, check C). In the last, a treated subject is at risk at a
  # control event, and the brentq root of the score is 0.8216793853198513.
  inf_sets <- list(
    list(1:6, c(1, 1, 1, 0, 0, 0), c(0, 0, 0, 1, 1, 1)),  # O_T = 0 = L
    list(1:4, c(1, 0, 1, 1), c(0, 0, 1, 1)),              # L = O_T = 2
    list(1:4, c(1, 1, 1, 0), c(1, 1, 0, 0))               # O_T = H = 2
  )
  for (d in inf_sets) {
    x <- coxph_fast(d[[1]], d[[2]], d[[3]], control = 0)
    expect_true(all(is.na(x)))
    expect_identical(attr(x, "converged"), NA)
  }
  fin <- coxph_fast(1:5, c(1, 1, 1, 0, 1), c(1, 1, 0, 1, 0), control = 0)
  expect_true(attr(fin, "converged"))
  expect_lt(abs(unname(fin["coef"]) - 0.8216793853198513), 1e-6)
  expect_equal(unname(fin["coef"]), 0.8216793837613253, tolerance = 1e-10)
})

test_that("coxph_fast converges when steps must be halved", {
  # Five treated subjects (n = 60), all with early events. Without the
  # bracket of the root kept after the first halved step, the iteration does
  # not converge in 20 iterations; with it, it converges in 8 with one halved
  # step (Planning/python_checks/cox_exact.py, check D, which also gives the
  # brentq root 4.901223249152088).
  tt <- c(0.6275, 0.6532, 0.2890, 0.1301, 0.0050, 0.2375, 0.0127, 0.0420,
          0.3368, 0.0600, 0.0034, 0.0951, 0.0783, 0.6512, 0.9587, 0.1691,
          0.4081, 0.7416, 0.2145, 0.3499, 0.1962, 0.9845, 0.1353, 0.3724,
          0.1657, 0.0214, 0.4140, 0.0436, 0.2786, 0.1745, 0.4734, 0.0043,
          0.8343, 0.4718, 0.7199, 0.2009, 0.0010, 0.0071, 0.0241, 0.0446,
          0.0235, 0.1780, 0.8053, 0.0619, 0.3778, 0.4031, 0.6465, 0.7278,
          0.1549, 0.2192, 0.1909, 0.8311, 0.0933, 0.0293, 0.3511, 0.8700,
          0.2604, 1.1072, 0.0729, 0.4103)
  ee <- c(0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0,
          1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1,
          0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0)
  xx <- c(0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  fit <- coxph_fast(tt, ee, xx, control = 0, ties = "breslow")
  expect_true(attr(fit, "converged"))
  expect_identical(attr(fit, "iter"), 8L)
  expect_equal(unname(fit["coef"]), 4.901223249152084, tolerance = 1e-10)
  expect_equal(unname(fit["se(coef)"]), 1.123500228719776, tolerance = 1e-10)
  expect_lt(abs(unname(fit["coef"]) - 4.901223249152088), 1e-6)
})

test_that("coxph_fast returns NA with a warning when the iteration fails", {
  # The core is replaced by a stub that reports no convergence
  local_mocked_bindings(cox_core = function(...) c(0.1, 0.2, -9, -8, 20, 1))
  expect_warning(
    x <- coxph_fast(1:6, c(1, 1, 1, 1, 0, 1), c(0, 1, 0, 1, 0, 1),
                    control = 0),
    "did not converge in 20 iterations")
  expect_true(all(is.na(x)))
  expect_false(attr(x, "converged"))
  expect_identical(attr(x, "iter"), 20L)
})

test_that("coxph_fast records the fit and validates ties", {
  set.seed(9)
  n  <- 80
  tt <- rexp(n, 0.1)                    # distinct times: no ties
  ee <- rbinom(n, 1, 0.8)
  gg <- rep(1:2, each = n / 2)
  ef <- coxph_fast(tt, ee, gg, control = 1)
  br <- coxph_fast(tt, ee, gg, control = 1, ties = "breslow")
  expect_identical(attr(ef, "ties"), "efron")
  expect_identical(attr(br, "ties"), "breslow")
  expect_true(is.integer(attr(ef, "iter")) && attr(ef, "iter") >= 1L)
  expect_true(attr(ef, "converged"))
  # Without tied event times the two handlings give the same result
  expect_identical(as.numeric(ef), as.numeric(br))
  expect_error(coxph_fast(tt, ee, gg, control = 1, ties = "exact"),
               "should be one of")
})

test_that("coxph_fast returns all-NA when no events", {
  time  <- c(1, 2, 3, 4)
  event <- c(0, 0, 0, 0)
  group <- c(1, 1, 2, 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_true(all(is.na(res)))
})

test_that("coxph_fast stops when input lengths differ", {
  expect_error(
    coxph_fast(c(1, 2, 3), c(1, 0), c(1, 1, 2), control = 1),
    "same length"
  )
})

test_that("coxph_fast validates event coding, group levels, and conf.level", {
  tt <- c(1, 2, 3, 4, 5, 6)
  gg <- c(1, 1, 1, 2, 2, 2)
  expect_error(coxph_fast(tt, c(1, 2, 1, 2, 1, 2), gg, control = 1),
               "coded as 0")
  expect_error(coxph_fast(tt, rep(1, 6), c(1, 1, 2, 2, 3, 3), control = 1),
               "two distinct")
  expect_error(coxph_fast(tt, rep(1, 6), gg, control = 3), "control")
  expect_error(coxph_fast(tt, rep(1, 6), gg, control = 1, conf.level = 95),
               "conf.level")
})

# ------------------------------------------------------------------ #
#  Stratified estimator
# ------------------------------------------------------------------ #

test_that("cox_core reproduces a stratified hand calculation", {
  # Stratum 1: (t, group, event) = (1, C, 1), (2, T, 1), (3, C, 0), (4, T, 1)
  # Stratum 2: (1, T, 1), (2, C, 1), (3, T, 1)
  # Risk sets restart in each stratum. At beta = 0, the sizes of the risk
  # sets give l(0) = -log(4 * 3 * 3 * 2) = -log(72), U(0) = O_T - E_T =
  # 4 - 13/3 = -1/3, and I(0) = 1/4 + 2/9 + 2/9 + 1/4 = 17/18, so the first
  # Newton-Raphson iterate is U(0) / I(0) = -6/17. With r = exp(beta) the
  # score is U = 2 - 2r / (r + 1) - 4r / (2r + 1), which vanishes at
  # r = 1 / sqrt(2), i.e. beta = -log(2) / 2 (Planning/python_checks/
  # cox_exact.py, check B).
  tt <- c(1, 2, 3, 4, 1, 2, 3)
  ee <- c(1L, 1L, 0L, 1L, 1L, 1L, 1L)
  jj <- c(0L, 1L, 0L, 1L, 1L, 0L, 1L)
  ss <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L)
  one <- cox_core(tt, ee, jj, ss, efron = FALSE, iter_max = 1)
  expect_equal(one[3], -log(72), tolerance = 1e-14)
  expect_equal(one[1], -6 / 17, tolerance = 1e-14)
  expect_identical(one[c(5, 6)], c(1, 1))
  full <- cox_core(tt, ee, jj, ss, efron = FALSE)
  expect_identical(full[6], 0)
  expect_lt(abs(full[1] + log(2) / 2), 1e-12)
  fast <- coxph_fast(tt, ee, jj, control = 0, strata = ss)
  expect_identical(unname(fast["coef"]), full[1])
})

test_that("a single stratum gives the unstratified estimate", {
  set.seed(401)
  n  <- 200
  tt <- rexp(n, 0.1)
  ee <- rbinom(n, 1, 0.8)
  gg <- rep(1:2, each = n / 2)
  un <- coxph_fast(tt, ee, gg, control = 1)
  st <- coxph_fast(tt, ee, gg, control = 1, strata = rep("all", n))
  expect_equal(as.numeric(st), as.numeric(un), tolerance = 1e-12)
  expect_equal(attr(st, "strata"), 1L)
})

test_that("stratified coxph_fast agrees with survival::coxph(strata())", {
  skip_if_not_installed("survival")
  set.seed(402)
  n_s  <- 150
  base <- c(0.05, 0.10, 0.20, 0.40)
  s    <- rep(seq_along(base), each = n_s)
  g    <- rep(0:1, times = length(base) * n_s / 2)
  lam  <- base[s] * ifelse(g == 1, 0.7, 1)
  t0   <- rexp(length(s), lam)
  cc   <- rexp(length(s), 0.05)
  tt   <- round(pmin(t0, cc), 1)       # rounding creates tied times
  ee   <- as.integer(t0 <= cc)
  ctrl <- survival::coxph.control(timefix = FALSE)
  for (tie in c("efron", "breslow")) {
    fast <- coxph_fast(tt, ee, g, control = 0, strata = s, ties = tie)
    ref  <- survival::coxph(survival::Surv(tt, ee) ~ g + survival::strata(s),
                            ties = tie, control = ctrl)
    expect_equal(unname(fast["coef"]), unname(stats::coef(ref)),
                 tolerance = 1e-10, label = tie)
    expect_equal(unname(fast["se(coef)"]), unname(sqrt(ref$var[1, 1])),
                 tolerance = 1e-10, label = tie)
    expect_equal(attr(fast, "strata"), 4L)
  }
})

test_that("stratified coxph_fast does not depend on the row order", {
  set.seed(403)
  n  <- 240
  s  <- sample(c("a", "b", "c"), n, replace = TRUE)
  tt <- rexp(n, ifelse(s == "a", 0.05, 0.2))
  ee <- rbinom(n, 1, 0.8)
  gg <- sample(1:2, n, replace = TRUE)
  ref <- coxph_fast(tt, ee, gg, control = 1, strata = s)
  shuf <- sample(n)
  expect_equal(as.numeric(coxph_fast(tt[shuf], ee[shuf], gg[shuf],
                                     control = 1, strata = s[shuf])),
               as.numeric(ref), tolerance = 1e-12)
  ord <- order(s, tt)
  expect_equal(as.numeric(coxph_fast(tt[ord], ee[ord], gg[ord], control = 1,
                                     strata = s[ord], presorted = TRUE)),
               as.numeric(ref), tolerance = 1e-12)
  expect_output(print(ref), "Stratified Cox proportional hazards model")
})

test_that("stratified coxph_fast validates strata", {
  tt <- c(1, 2, 3, 4, 5, 6)
  gg <- c(1, 1, 1, 2, 2, 2)
  ee <- rep(1, 6)
  expect_error(coxph_fast(tt, ee, gg, control = 1, strata = c(1, 2)),
               "same length")
  expect_error(coxph_fast(tt, ee, gg, control = 1,
                          strata = c(1, NA, 1, 2, 2, 2)), "missing")
  expect_error(coxph_fast(tt, ee, gg, control = 1, presorted = TRUE,
                          strata = c(1, 2, 1, 2, 1, 2)), "contiguous")
})

test_that("coxph_fast: presorted = TRUE checks the order", {
  time  <- c(3, 1, 2, 4, 6, 5)
  event <- c(1, 1, 0, 1, 1, 0)
  group <- c(0, 1, 0, 1, 0, 1)
  expect_error(coxph_fast(time, event, group, control = 0, presorted = TRUE), "presorted = FALSE")
  expect_error(coxph_fast(time, event, group, control = 0, presorted = TRUE, strata = c(1, 1, 1, 2, 2, 2)), "presorted = FALSE")
  expect_error(coxph_fast(c(1, 3, 2, 4, 5, 6), event, group, control = 0, presorted = TRUE, strata = c(1, 2, 1, 2, 1, 2)), "presorted = FALSE")
})

test_that("coxph_fast: the formula interface gives the vector result", {
  set.seed(63)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1), ee = rbinom(n, 1, 0.8),
                  g = rep(c(1, 0), n / 2))
  d$s <- rep(c("x", "y", "z"), n / 3)
  ref <- coxph_fast(d$tt, d$ee, d$g, control = 0)
  # The control group defaults to the smallest value of g
  expect_identical(coxph_fast(Surv(tt, ee) ~ g, data = d), ref)
  expect_identical(coxph_fast(Surv(tt, ee) ~ g, d, side = 1),
                   coxph_fast(d$tt, d$ee, d$g, control = 0, side = 1))
  expect_identical(coxph_fast(Surv(tt, ee) ~ g, d, control = 1),
                   coxph_fast(d$tt, d$ee, d$g, control = 1))
  # An explicit NULL control also gives the first level
  expect_identical(coxph_fast(Surv(tt, ee) ~ g, d, control = NULL), ref)
  expect_identical(class(ref), c("coxph_fast", "fastsurv"))
  expect_identical(coxph_fast(Surv(tt, ee) ~ g + strata(s), data = d),
                   coxph_fast(d$tt, d$ee, d$g, control = 0, strata = d$s))
  expect_error(coxph_fast(Surv(tt, ee) ~ g, d, 0), "right-hand side")
  expect_error(coxph_fast(d$tt, d$ee, d$g, control = 0, data = d),
               "only with a formula")
})

test_that("coxph_fast: the formula interface passes every argument", {
  set.seed(60)
  n <- 120
  d <- data.frame(tt = round(rexp(n, 0.1), 1) + 0.1, ee = rbinom(n, 1, 0.8),
                  g = rep(c("b", "a"), n / 2), s = rep(c("x", "y", "z"), n / 3))
  d1 <- d[order(d$tt), ]
  d2 <- d[order(d$s, d$tt), ]
  expect_identical(
    coxph_fast(Surv(tt, ee) ~ g + strata(s), d2, control = "b", side = 1,
               conf.level = 0.9, presorted = TRUE, ties = "breslow"),
    coxph_fast(d2$tt, d2$ee, d2$g, "b", side = 1, conf.level = 0.9,
               presorted = TRUE, strata = d2$s, ties = "breslow"))
  # The control group defaults to the first value in sorted order
  expect_identical(coxph_fast(Surv(tt, ee) ~ g, d),
                   coxph_fast(d$tt, d$ee, d$g, control = "a"))
})

test_that("coxph_fast: a data frame can be piped into the formula interface", {
  set.seed(703)
  d <- data.frame(tt = round(rexp(120, 0.1), 1) + 0.1,
                  ee = rbinom(120, 1, 0.8), g = rep(c(1, 0), 60))
  fo <- Surv(tt, ee) ~ g
  expect_identical(d |> coxph_fast(fo),
                   coxph_fast(fo, data = d))
  expect_error(d |> coxph_fast(fo, data = d), "given twice")
})
