test_that("simdata_fast returns a data.frame with correct columns (one group)", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 20,
    a.time   = c(0, 12),
    a.rate   = 20 / 12,
    e.median = 18,
    seed     = 1
  )

  expect_s3_class(df, "data.frame")
  expect_named(df, c("sim", "group", "accrual_time", "surv_time",
                     "dropout_time", "tte", "event", "calendar_time"))
})

test_that("simdata_fast returns correct number of rows (one group)", {
  nsim <- 10
  n    <- 20
  df   <- simdata_fast(
    nsim     = nsim,
    n        = n,
    a.time   = c(0, 12),
    a.rate   = n / 12,
    e.median = 18,
    seed     = 2
  )

  expect_equal(nrow(df), nsim * n)
})

test_that("simdata_fast returns correct number of rows (two groups)", {
  nsim <- 10
  n    <- c(50, 50)
  df   <- simdata_fast(
    nsim     = nsim,
    n        = n,
    a.time   = c(0, 12),
    a.rate   = sum(n) / 12,
    e.median = list(18, 24),
    seed     = 3
  )

  expect_equal(nrow(df), nsim * sum(n))
  expect_equal(sort(unique(df$group)), c(1L, 2L))
})

test_that("simdata_fast accrual_time is within a.time range", {
  df <- simdata_fast(
    nsim     = 20,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    seed     = 4
  )

  expect_gte(min(df$accrual_time), 0)
  expect_lte(max(df$accrual_time), 12)
})

test_that("simdata_fast tte = pmin(surv_time, dropout_time)", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    d.hazard = 0.02,
    seed     = 5
  )

  expect_equal(df$tte, pmin(df$surv_time, df$dropout_time))
})

test_that("simdata_fast event = 1 iff surv_time <= dropout_time", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    d.hazard = 0.02,
    seed     = 6
  )

  expected_event <- as.integer(df$surv_time <= df$dropout_time)
  expect_equal(df$event, expected_event)
})

test_that("simdata_fast calendar_time = accrual_time + tte", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    seed     = 7
  )

  expect_equal(df$calendar_time, df$accrual_time + df$tte)
})

test_that("simdata_fast dropout_time = Inf when no dropout specified", {
  df <- simdata_fast(
    nsim     = 10,
    n        = 50,
    a.time   = c(0, 12),
    a.rate   = 50 / 12,
    e.median = 18,
    seed     = 8
  )

  expect_true(all(is.infinite(df$dropout_time)))
  expect_true(all(df$event == 1L))
})

test_that("simdata_fast is reproducible with seed", {
  df1 <- simdata_fast(
    nsim     = 5,
    n        = 20,
    a.time   = c(0, 12),
    a.rate   = 20 / 12,
    e.median = 18,
    seed     = 42
  )
  df2 <- simdata_fast(
    nsim     = 5,
    n        = 20,
    a.time   = c(0, 12),
    a.rate   = 20 / 12,
    e.median = 18,
    seed     = 42
  )

  expect_equal(df1, df2)
})

test_that("simdata_fast piecewise exponential: surv_time > 0", {
  df <- simdata_fast(
    nsim     = 10,
    n        = c(50, 50),
    a.time   = c(0, 12),
    a.rate   = 100 / 12,
    e.hazard = list(c(0.08, 0.08), c(0.08, 0.04)),
    e.time   = c(0, 6, Inf),
    seed     = 9
  )

  expect_true(all(df$surv_time > 0))
})

test_that("simdata_fast e.hazard and e.median are mutually exclusive", {
  expect_error(
    simdata_fast(
      nsim     = 5,
      n        = 20,
      a.time   = c(0, 12),
      a.rate   = 20 / 12,
      e.hazard = 0.05,
      e.median = 18,
      seed     = 10
    ),
    "exactly one"
  )
})

test_that("simdata_fast total n with alloc splits correctly", {
  nsim  <- 5
  n_tot <- 100
  df    <- simdata_fast(
    nsim     = nsim,
    n        = n_tot,
    alloc    = c(1, 1),
    a.time   = c(0, 12),
    a.rate   = n_tot / 12,
    e.median = list(18, 24),
    seed     = 11
  )

  n_per_sim <- table(df$sim[df$sim == 1 & df$group == 1])
  expect_equal(nrow(df), nsim * n_tot)
})

# ------------------------------------------------------------------ #
#  Subgroups: single factor
# ------------------------------------------------------------------ #

test_that("simdata_fast single-factor prevalence inserts a 'subgroup' column", {
  df <- simdata_fast(
    nsim       = 10,
    n          = 100,
    a.time     = c(0, 12),
    a.rate     = 100 / 12,
    e.hazard   = 0.05,
    prevalence = c(0.5, 0.3, 0.2),
    seed       = 12
  )

  expect_named(df, c("sim", "group", "subgroup", "accrual_time", "surv_time",
                     "dropout_time", "tte", "event", "calendar_time"))
  expect_true(is.integer(df$subgroup))
  expect_true(all(df$subgroup %in% 1:3))
})

test_that("simdata_fast prevalence = NULL matches degenerate single subgroup", {
  args <- list(nsim = 10, n = c(50, 60), a.time = c(0, 12),
               a.rate = 110 / 12, e.median = list(18, 24),
               d.hazard = list(0.01, 0.02), seed = 13)

  df_null <- do.call(simdata_fast, args)
  df_one  <- do.call(simdata_fast, c(args, list(prevalence = c(1))))

  # The degenerate K = 1 path takes the same numeric route (no subgroup
  # column, no extra random draws), so the simulated values coincide.
  expect_equal(as.list(df_null), as.list(df_one))
})

test_that("simdata_fast single-factor prevalence recovers proportions", {
  df <- simdata_fast(
    nsim       = 1,
    n          = 50000,
    a.time     = c(0, 12),
    a.rate     = 50000 / 12,
    e.hazard   = 0.05,
    prevalence = c(0.5, 0.3, 0.2),
    seed       = 14
  )

  prop <- as.numeric(table(factor(df$subgroup, levels = 1:3))) / nrow(df)
  expect_equal(prop, c(0.5, 0.3, 0.2), tolerance = 0.02)
})

test_that("simdata_fast per-subgroup hazard recovers mean survival", {
  haz <- list(0.10, 0.05)
  df  <- simdata_fast(
    nsim       = 1,
    n          = 50000,
    a.time     = c(0, 12),
    a.rate     = 50000 / 12,
    e.hazard   = haz,
    prevalence = c(0.5, 0.5),
    seed       = 15
  )

  m1 <- mean(df$surv_time[df$subgroup == 1L])
  m2 <- mean(df$surv_time[df$subgroup == 2L])
  expect_equal(m1, 1 / haz[[1L]], tolerance = 0.05)
  expect_equal(m2, 1 / haz[[2L]], tolerance = 0.05)
})

# ------------------------------------------------------------------ #
#  Subgroups: multiple factors
# ------------------------------------------------------------------ #

test_that("simdata_fast independent factors insert subgroup1/subgroup2", {
  df <- simdata_fast(
    nsim       = 10,
    n          = 200,
    a.time     = c(0, 12),
    a.rate     = 200 / 12,
    e.hazard   = 0.05,
    prevalence = list(c(0.5, 0.5), c(0.6, 0.4)),
    seed       = 16
  )

  expect_named(df, c("sim", "group", "subgroup1", "subgroup2", "accrual_time",
                     "surv_time", "dropout_time", "tte", "event",
                     "calendar_time"))
  expect_true(all(df$subgroup1 %in% 1:2))
  expect_true(all(df$subgroup2 %in% 1:2))
})

test_that("simdata_fast independent factors are independent", {
  m1 <- c(0.7, 0.3)
  m2 <- c(0.5, 0.3, 0.2)
  df <- simdata_fast(
    nsim       = 1,
    n          = 100000,
    a.time     = c(0, 12),
    a.rate     = 100000 / 12,
    e.hazard   = 0.05,
    prevalence = list(m1, m2),
    seed       = 17
  )

  joint_emp  <- table(factor(df$subgroup1, levels = 1:2),
                      factor(df$subgroup2, levels = 1:3)) / nrow(df)
  joint_theo <- outer(m1, m2)
  expect_lt(max(abs(joint_emp - joint_theo)), 0.01)
})

test_that("simdata_fast joint-distribution array captures dependence", {
  jd <- array(c(0.40, 0.10, 0.15, 0.35), dim = c(2, 2))
  df <- simdata_fast(
    nsim       = 1,
    n          = 100000,
    a.time     = c(0, 12),
    a.rate     = 100000 / 12,
    e.hazard   = 0.05,
    prevalence = jd,
    seed       = 18
  )

  joint_emp <- table(factor(df$subgroup1, levels = 1:2),
                     factor(df$subgroup2, levels = 1:2)) / nrow(df)
  expect_lt(max(abs(as.numeric(joint_emp) - as.numeric(jd))), 0.01)

  # The empirical joint should not match the product of marginals.
  m1 <- rowSums(jd); m2 <- colSums(jd)
  expect_gt(max(abs(as.numeric(joint_emp) - as.numeric(outer(m1, m2)))), 0.05)
})

test_that("simdata_fast per-cell hazard uses column-major cell order", {
  # Two 2-level factors -> cells (f1,f2): 1=(1,1), 2=(2,1), 3=(1,2), 4=(2,2)
  haz <- list(0.10, 0.05, 0.04, 0.02)
  df  <- simdata_fast(
    nsim       = 1,
    n          = 80000,
    a.time     = c(0, 12),
    a.rate     = 80000 / 12,
    e.hazard   = haz,
    prevalence = list(c(0.5, 0.5), c(0.5, 0.5)),
    seed       = 19
  )

  cell_of <- function(f1, f2) (f2 - 1L) * 2L + f1
  for (f1 in 1:2) for (f2 in 1:2) {
    cidx <- cell_of(f1, f2)
    m    <- mean(df$surv_time[df$subgroup1 == f1 & df$subgroup2 == f2])
    expect_equal(m, 1 / haz[[cidx]], tolerance = 0.05)
  }
})

test_that("simdata_fast group-specific prevalence recovers per-group composition", {
  prev <- list(control = c(0.7, 0.3), treatment = c(0.4, 0.6))
  df <- simdata_fast(
    nsim       = 1,
    n          = c(50000, 50000),
    a.time     = c(0, 12),
    a.rate     = 100000 / 12,
    e.hazard   = list(0.05, 0.05),
    prevalence = prev,
    seed       = 20
  )

  g1 <- df[df$group == 1L, ]
  g2 <- df[df$group == 2L, ]
  p1 <- as.numeric(table(factor(g1$subgroup, levels = 1:2))) / nrow(g1)
  p2 <- as.numeric(table(factor(g2$subgroup, levels = 1:2))) / nrow(g2)
  expect_equal(p1, prev$control,   tolerance = 0.02)
  expect_equal(p2, prev$treatment, tolerance = 0.02)
})

test_that("simdata_fast group-specific prevalence errors on mismatched levels", {
  expect_error(
    simdata_fast(
      nsim       = 5,
      n          = c(100, 100),
      a.time     = c(0, 12),
      a.rate     = 200 / 12,
      e.hazard   = list(0.05, 0.05),
      prevalence = list(control = c(0.5, 0.5),
                        treatment = c(0.3, 0.3, 0.4)),
      seed       = 21
    ),
    "same"
  )
})

test_that("simdata_fast errors on negative prevalence", {
  expect_error(
    simdata_fast(
      nsim       = 5,
      n          = 100,
      a.time     = c(0, 12),
      a.rate     = 100 / 12,
      e.hazard   = 0.05,
      prevalence = c(0.5, -0.1, 0.6),
      seed       = 22
    ),
    "positive"
  )
})

test_that("simdata_fast errors on per-cell hazard list of wrong length", {
  expect_error(
    simdata_fast(
      nsim       = 5,
      n          = 200,
      a.time     = c(0, 12),
      a.rate     = 200 / 12,
      e.hazard   = list(0.1, 0.05, 0.04),  # 3 specs for 4 cells
      prevalence = list(c(0.5, 0.5), c(0.5, 0.5)),
      seed       = 23
    ),
    "per cell"
  )
})

test_that("simdata_fast normalizes unnormalized prevalence", {
  raw  <- c(2, 1, 1)
  norm <- raw / sum(raw)
  df <- simdata_fast(
    nsim       = 1,
    n          = 50000,
    a.time     = c(0, 12),
    a.rate     = 50000 / 12,
    e.hazard   = 0.05,
    prevalence = raw,
    seed       = 24
  )

  prop <- as.numeric(table(factor(df$subgroup, levels = 1:3))) / nrow(df)
  expect_equal(prop, norm, tolerance = 0.02)
})

# ------------------------------------------------------------------ #
#  Fixed (deterministic) allocation
# ------------------------------------------------------------------ #

test_that("simdata_fast fixed.alloc n=99 c(.5,.5) gives 50/49 every sim", {
  df <- simdata_fast(
    nsim        = 10,
    n           = 99,
    a.time      = c(0, 12),
    a.rate      = 99 / 12,
    e.hazard    = 0.05,
    prevalence  = c(0.5, 0.5),
    fixed.alloc = TRUE,
    seed        = 25
  )

  counts <- tapply(df$subgroup, df$sim,
                   function(x) as.integer(table(factor(x, levels = 1:2))))
  m <- do.call(rbind, counts)
  expect_true(all(m[, 1] == 50L))
  expect_true(all(m[, 2] == 49L))
})

test_that("simdata_fast fixed.alloc n=100 thirds gives 34/33/33 every sim", {
  df <- simdata_fast(
    nsim        = 10,
    n           = 100,
    a.time      = c(0, 12),
    a.rate      = 100 / 12,
    e.hazard    = 0.05,
    prevalence  = c(1/3, 1/3, 1/3),
    fixed.alloc = TRUE,
    seed        = 26
  )

  counts <- tapply(df$subgroup, df$sim,
                   function(x) as.integer(table(factor(x, levels = 1:3))))
  m <- do.call(rbind, counts)
  expect_true(all(m[, 1] == 34L))
  expect_true(all(m[, 2] == 33L))
  expect_true(all(m[, 3] == 33L))
})

test_that("simdata_fast fixed.alloc cell counts always sum to n", {
  df <- simdata_fast(
    nsim        = 8,
    n           = 97,
    a.time      = c(0, 12),
    a.rate      = 97 / 12,
    e.hazard    = 0.05,
    prevalence  = c(0.4, 0.35, 0.25),
    fixed.alloc = TRUE,
    seed        = 27
  )

  counts <- tapply(df$subgroup, df$sim, length)
  expect_true(all(counts == 97L))
})

test_that("simdata_fast fixed.alloc labels are deterministic across seeds", {
  args <- list(nsim = 5, n = 97, a.time = c(0, 12), a.rate = 97 / 12,
               e.hazard = 0.05, prevalence = c(0.4, 0.35, 0.25),
               fixed.alloc = TRUE)
  da <- do.call(simdata_fast, c(args, list(seed = 28)))
  db <- do.call(simdata_fast, c(args, list(seed = 999)))

  # Cell labels do not depend on the RNG seed, but survival times do.
  expect_equal(da$subgroup, db$subgroup)
  expect_false(isTRUE(all.equal(da$surv_time, db$surv_time)))
})

test_that("simdata_fast fixed.alloc multi-factor cell counts are exact", {
  # cell_prob column-major = (.5*.6, .5*.6, .5*.4, .5*.4) = (.30,.30,.20,.20)
  # n = 10 -> (3, 3, 2, 2)
  df <- simdata_fast(
    nsim        = 5,
    n           = 10,
    a.time      = c(0, 12),
    a.rate      = 10 / 12,
    e.hazard    = 0.05,
    prevalence  = list(c(0.5, 0.5), c(0.6, 0.4)),
    fixed.alloc = TRUE,
    seed        = 29
  )

  for (s in 1:5) {
    d    <- df[df$sim == s, ]
    cidx <- (d$subgroup2 - 1L) * 2L + d$subgroup1
    cnt  <- as.integer(table(factor(cidx, levels = 1:4)))
    expect_equal(cnt, c(3L, 3L, 2L, 2L))
  }
})
