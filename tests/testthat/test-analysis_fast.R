# Reference administrative censoring at a calendar cutoff, optionally
# restricted to a subgroup mask, returned unsorted. Mirrors the logic in
# analysis_cut_core so the R-level plumbing can be checked independently.
manual_cut <- function(acc, tte, ev, j, cutoff, mask = NULL) {
  enrolled <- acc <= cutoff
  if (!is.null(mask)) enrolled <- enrolled & mask
  a    <- acc[enrolled]
  full <- tte[enrolled]
  e0   <- ev[enrolled]
  jj   <- j[enrolled]
  before <- (a + full) <= cutoff
  list(time = ifelse(before, full, cutoff - a),
       event = as.integer(ifelse(before, e0, 0L)),
       j = as.integer(jj),
       n = sum(enrolled),
       n_event = sum(ifelse(before, e0, 0L)))
}

make_data <- function(seed = 101, nsim = 30) {
  simdata_fast(
    nsim       = nsim,
    n          = c(120, 120),
    a.time     = c(0, 12),
    a.rate     = 240 / 12,
    e.hazard   = list(list(0.05, 0.035), 0.030),
    prevalence = c(0.5, 0.5),
    seed       = seed
  )
}

# ------------------------------------------------------------------ #
#  Structure and input validation
# ------------------------------------------------------------------ #

test_that("analysis_fast returns one row per (sim, look) without subgroups", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, time.looks = c(18, 30))
  expect_s3_class(res, "data.frame")
  expect_false("population" %in% names(res))
  expect_equal(nrow(res), 30L * 2L)
})

test_that("analysis_fast errors when both look types supplied", {
  dat <- make_data()
  expect_error(
    analysis_fast(dat, control = 1, event.looks = 80, time.looks = 24),
    "exactly one"
  )
})

test_that("analysis_fast errors when neither look type supplied", {
  dat <- make_data()
  expect_error(analysis_fast(dat, control = 1), "exactly one")
})

test_that("analysis_fast errors on rmst without tau and km without t.eval", {
  dat <- make_data()
  expect_error(
    analysis_fast(dat, control = 1, time.looks = 24, stat = "rmst"),
    "tau"
  )
  expect_error(
    analysis_fast(dat, control = 1, time.looks = 24, stat = "km"),
    "t.eval"
  )
})

test_that("analysis_fast errors on invalid side", {
  dat <- make_data()
  expect_error(
    analysis_fast(dat, control = 1, time.looks = 24, side = 3),
    "side"
  )
})

# ------------------------------------------------------------------ #
#  Calendar look: agreement with a manual administrative censoring
# ------------------------------------------------------------------ #

test_that("analysis_fast calendar look counts match manual censoring", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, time.looks = c(18, 30),
                       stat = "logrank")
  for (s in c(1, 2, 5)) {
    rows <- which(dat$sim == s)
    acc  <- dat$accrual_time[rows]; tt <- dat$tte[rows]
    ev   <- as.integer(dat$event[rows]); jj <- as.integer(dat$group[rows] != 1)
    for (l in 1:2) {
      mc  <- manual_cut(acc, tt, ev, jj, c(18, 30)[l])
      row <- res[res$sim == s & res$look == l, ]
      expect_equal(row$n.enrolled, mc$n)
      expect_equal(row$n.event, mc$n_event)
    }
  }
})

test_that("analysis_fast log-rank/Cox/RMST match direct calls on cut data", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, time.looks = 30,
                       stat = c("logrank", "coxph", "rmst"), tau = 15)
  for (s in c(1, 3, 8)) {
    rows <- which(dat$sim == s)
    acc  <- dat$accrual_time[rows]; tt <- dat$tte[rows]
    ev   <- as.integer(dat$event[rows]); jj <- as.integer(dat$group[rows] != 1)
    mc   <- manual_cut(acc, tt, ev, jj, 30)
    row  <- res[res$sim == s, ]

    z_ref <- as.numeric(survdiff_fast(mc$time, mc$event, mc$j,
                                      control = 0L, side = 1L))
    expect_equal(row$logrank.z, z_ref, tolerance = 1e-8)

    cx_ref <- coxph_fast(mc$time, mc$event, mc$j, control = 0L)
    expect_equal(row$cox.coef, unname(cx_ref[1L]), tolerance = 1e-8)
    expect_equal(row$cox.hr,   unname(cx_ref[2L]), tolerance = 1e-8)

    rm_ref <- rmst_fast(mc$time, mc$event, group = mc$j, control = 0L, tau = 15)
    expect_equal(row$rmst.diff, unname(rm_ref["diff"]), tolerance = 1e-8)
  }
})

# ------------------------------------------------------------------ #
#  Event-driven look
# ------------------------------------------------------------------ #

test_that("analysis_fast event look stops at the d-th event time", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, event.looks = c(50, 120),
                       stat = "logrank")
  for (s in c(1, 2, 5)) {
    rows   <- which(dat$sim == s)
    acc    <- dat$accrual_time[rows]; tt <- dat$tte[rows]
    ev     <- as.integer(dat$event[rows])
    cal_ev <- (acc + tt)[ev == 1L]
    for (l in 1:2) {
      d   <- c(50, 120)[l]
      row <- res[res$sim == s & res$look == l, ]
      if (length(cal_ev) >= d) {
        expect_equal(row$cutoff, sort(cal_ev)[d], tolerance = 1e-9)
        expect_equal(row$n.event, d)
        expect_true(row$reached)
      }
    }
  }
})

test_that("analysis_fast unreached event target falls back to full data", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, event.looks = 100000L,
                       stat = "logrank")
  expect_true(all(!res$reached))
  expect_true(all(is.na(res$cutoff)))

  for (s in c(1, 3)) {
    rows <- which(dat$sim == s)
    tt   <- dat$tte[rows]; ev <- as.integer(dat$event[rows])
    jj   <- as.integer(dat$group[rows] != 1)
    row  <- res[res$sim == s, ]
    z_full <- as.numeric(survdiff_fast(tt, ev, jj, control = 0L, side = 1L))
    expect_equal(row$logrank.z, z_full, tolerance = 1e-8)
    expect_equal(row$n.enrolled, length(tt))
    expect_equal(row$n.event, sum(ev))
  }
})

# ------------------------------------------------------------------ #
#  side and p-values
# ------------------------------------------------------------------ #

test_that("analysis_fast two-sided p-values equal 2 * pnorm(-|z|)", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, time.looks = 30,
                       stat = c("logrank", "coxph", "rmst"),
                       tau = 15, side = 2)
  ok  <- !is.na(res$logrank.z)
  expect_equal(res$logrank.p[ok], 2 * pnorm(-abs(res$logrank.z[ok])),
               tolerance = 1e-10)
  expect_equal(res$cox.p[ok], 2 * pnorm(-abs(res$cox.z[ok])),
               tolerance = 1e-10)
  expect_equal(res$rmst.p[ok], 2 * pnorm(-abs(res$rmst.z[ok])),
               tolerance = 1e-10)
})

test_that("analysis_fast one-sided p-values use each test's benefit direction", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, time.looks = 30,
                       stat = c("logrank", "coxph", "rmst"),
                       tau = 15, side = 1)
  ok <- !is.na(res$logrank.z)
  # log-rank and Cox: benefit is a negative Z -> lower tail pnorm(z)
  expect_equal(res$logrank.p[ok], pnorm(res$logrank.z[ok]), tolerance = 1e-10)
  expect_equal(res$cox.p[ok], pnorm(res$cox.z[ok]), tolerance = 1e-10)
  # RMST: benefit is a positive Z -> upper tail pnorm(-z)
  expect_equal(res$rmst.p[ok], pnorm(-res$rmst.z[ok]), tolerance = 1e-10)
})

test_that("analysis_fast Cox z matches coef/se and shares the coef sign", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, time.looks = 30, stat = "coxph")
  ok  <- !is.na(res$cox.z)
  expect_equal(res$cox.z[ok], res$cox.coef[ok] / res$cox.se[ok],
               tolerance = 1e-10)
  # Natural sign: log HR and its Wald z share the same sign.
  expect_true(all(sign(res$cox.z[ok]) == sign(res$cox.coef[ok])))
})

test_that("analysis_fast log-rank and Cox z point the same direction", {
  # Both are tests on the log hazard ratio, so their signed Z should agree
  # in sign on most simulations (allowing rare boundary disagreements).
  dat <- make_data(seed = 202, nsim = 100)
  res <- analysis_fast(dat, control = 1, time.looks = 30,
                       stat = c("logrank", "coxph"))
  ok  <- !is.na(res$logrank.z) & !is.na(res$cox.z)
  agree <- mean(sign(res$logrank.z[ok]) == sign(res$cox.z[ok]))
  expect_gt(agree, 0.95)
})

# ------------------------------------------------------------------ #
#  Subgroup output (long format)
# ------------------------------------------------------------------ #

test_that("analysis_fast by.subgroup adds a population column and rows", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, time.looks = c(18, 30),
                       stat = "logrank", by.subgroup = TRUE)
  expect_true("population" %in% names(res))
  expect_setequal(unique(res$population),
                  c("overall", "subgroup_1", "subgroup_2"))
  expect_equal(nrow(res), 30L * 2L * 3L)
})

test_that("analysis_fast overall rows equal the whole-population analysis", {
  dat   <- make_data()
  res_s <- analysis_fast(dat, control = 1, time.looks = c(18, 30),
                         stat = c("logrank", "coxph"), by.subgroup = TRUE)
  res_p <- analysis_fast(dat, control = 1, time.looks = c(18, 30),
                         stat = c("logrank", "coxph"))
  ov <- res_s[res_s$population == "overall", ]
  for (s in c(1, 2, 5)) for (l in 1:2) {
    a <- ov[ov$sim == s & ov$look == l, ]
    b <- res_p[res_p$sim == s & res_p$look == l, ]
    expect_equal(a$logrank.z, b$logrank.z, tolerance = 1e-12)
    expect_equal(a$cox.coef, b$cox.coef, tolerance = 1e-12)
    expect_equal(a$n.event, b$n.event)
  }
})

test_that("analysis_fast subgroup counts partition the overall counts", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, time.looks = 30,
                       stat = "logrank", by.subgroup = TRUE)
  for (s in c(1, 2, 5)) {
    sub <- res[res$sim == s, ]
    ov_n <- sub$n.enrolled[sub$population == "overall"]
    ov_e <- sub$n.event[sub$population == "overall"]
    s_n  <- sum(sub$n.enrolled[sub$population != "overall"])
    s_e  <- sum(sub$n.event[sub$population != "overall"])
    expect_equal(s_n, ov_n)
    expect_equal(s_e, ov_e)
  }
})

test_that("analysis_fast subgroup rows match manual subgroup censoring", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, time.looks = 30,
                       stat = c("logrank", "coxph"), by.subgroup = TRUE)
  for (s in c(1, 3, 8)) {
    rows <- which(dat$sim == s)
    acc  <- dat$accrual_time[rows]; tt <- dat$tte[rows]
    ev   <- as.integer(dat$event[rows]); jj <- as.integer(dat$group[rows] != 1)
    sg   <- dat$subgroup[rows]
    for (lev in 1:2) {
      mc  <- manual_cut(acc, tt, ev, jj, 30, mask = sg == lev)
      row <- res[res$sim == s & res$population == paste0("subgroup_", lev), ]
      expect_equal(row$n.enrolled, mc$n)
      expect_equal(row$n.event, mc$n_event)
      if (mc$n_event > 0 && sum(mc$j == 0) > 0 && sum(mc$j == 1) > 0) {
        z_ref <- as.numeric(survdiff_fast(mc$time, mc$event, mc$j,
                                          control = 0L, side = 1L))
        expect_equal(row$logrank.z, z_ref, tolerance = 1e-8)
      }
    }
  }
})

test_that("analysis_fast cutoff and reached are shared across populations", {
  dat <- make_data()
  res <- analysis_fast(dat, control = 1, event.looks = 80,
                       stat = "logrank", by.subgroup = TRUE)
  for (s in c(1, 2, 5)) {
    sub <- res[res$sim == s, ]
    expect_equal(length(unique(sub$cutoff)), 1L)
    expect_equal(length(unique(sub$reached)), 1L)
  }
})

test_that("analysis_fast multi-factor subgroups label both factors", {
  dat <- simdata_fast(nsim = 20, n = c(100, 100), a.time = c(0, 12),
                      a.rate = 200 / 12, e.hazard = list(0.04, 0.03),
                      prevalence = list(c(0.5, 0.5), c(0.6, 0.4)), seed = 31)
  res <- analysis_fast(dat, control = 1, time.looks = 24,
                       stat = "logrank", by.subgroup = TRUE)
  expect_setequal(unique(res$population),
                  c("overall", "subgroup1_1", "subgroup1_2",
                    "subgroup2_1", "subgroup2_2"))
})

test_that("analysis_fast errors with by.subgroup when no subgroup columns", {
  dat <- simdata_fast(nsim = 5, n = c(50, 50), a.time = c(0, 12),
                      a.rate = 100 / 12, e.median = list(18, 24), seed = 32)
  expect_error(
    analysis_fast(dat, control = 1, time.looks = 24, by.subgroup = TRUE),
    "no subgroup"
  )
})

# ------------------------------------------------------------------ #
#  External numerical agreement
# ------------------------------------------------------------------ #

test_that("analysis_fast log-rank chi-square matches survival::survdiff", {
  skip_if_not_installed("survival")
  dat  <- make_data()
  res  <- analysis_fast(dat, control = 1, time.looks = 30, stat = "logrank")
  rows <- which(dat$sim == 1)
  acc  <- dat$accrual_time[rows]; tt <- dat$tte[rows]
  ev   <- as.integer(dat$event[rows]); jj <- as.integer(dat$group[rows] != 1)
  mc   <- manual_cut(acc, tt, ev, jj, 30)

  ref <- survival::survdiff(survival::Surv(mc$time, mc$event) ~ mc$j)
  row <- res[res$sim == 1, ]
  expect_equal(row$logrank.chisq, ref$chisq, tolerance = 1e-6)
})

test_that("analysis_fast Cox coef is close to survival::coxph (Breslow)", {
  skip_if_not_installed("survival")
  dat  <- make_data()
  res  <- analysis_fast(dat, control = 1, time.looks = 30, stat = "coxph")
  rows <- which(dat$sim == 1)
  acc  <- dat$accrual_time[rows]; tt <- dat$tte[rows]
  ev   <- as.integer(dat$event[rows]); jj <- as.integer(dat$group[rows] != 1)
  mc   <- manual_cut(acc, tt, ev, jj, 30)

  fit <- survival::coxph(survival::Surv(mc$time, mc$event) ~ mc$j,
                         ties = "breslow")
  row <- res[res$sim == 1, ]
  expect_equal(row$cox.coef, unname(stats::coef(fit)), tolerance = 5e-3)
})

test_that("analysis_fast RMST difference matches survRM2", {
  skip_if_not_installed("survRM2")
  dat  <- make_data()
  res  <- analysis_fast(dat, control = 1, time.looks = 30,
                        stat = "rmst", tau = 15)
  rows <- which(dat$sim == 1)
  acc  <- dat$accrual_time[rows]; tt <- dat$tte[rows]
  ev   <- as.integer(dat$event[rows]); jj <- as.integer(dat$group[rows] != 1)
  mc   <- manual_cut(acc, tt, ev, jj, 30)
  skip_if(any(mc$time <= 0))

  ref <- survRM2::rmst2(time = mc$time, status = mc$event, arm = mc$j, tau = 15)
  row <- res[res$sim == 1, ]
  expect_equal(row$rmst.diff, ref$unadjusted.result[1L, 1L], tolerance = 1e-6)
})
