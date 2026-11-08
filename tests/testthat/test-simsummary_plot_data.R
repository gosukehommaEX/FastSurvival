# Analysis output of six simulated trials at two looks. The expected values
# were computed independently from the rules in ?simsummary_fast and
# ?plot.simsummary_fast with exact fractions
# (Planning/python_checks/simsummary_plot_data.py of the development
# repository).
spd_data <- function() {
  df <- data.frame(
    sim       = rep(1:6, each = 2),
    look      = rep(1:2, 6),
    logrank.z = c(-2.9, -3.1, -1.0, -2.2, 0.5, -1.0,
                  -2.0, -2.5, -0.2, 0.3, -1.5, -2.1),
    cox.coef  = c(-0.4, -0.5, -0.1, -0.3, 0.3, 0.1,
                  -0.2, -0.3, 0.25, 0.2, -0.2, -0.3),
    n.event   = rep(c(50, 100), 6),
    cutoff    = rep(c(12, 24), 6)
  )
  df$logrank.p <- pnorm(df$logrank.z)
  df
}

test_that("simsummary_plot_data gives the stopping probabilities drawn", {
  df  <- spd_data()
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96),
                         fut.col = "cox.coef", futility = c(0.2, NA))

  pd <- simsummary_plot_data(res)
  expect_true(pd$has.futility)
  expect_null(pd$xcol)
  expect_identical(pd$looks, c("1", "2"))
  expect_length(pd$blocks, 1L)
  b <- pd$blocks[[1L]]
  expect_identical(b$label, "")
  expect_equal(b$x, c(1, 2))
  expect_equal(b$cum.efficacy, c(1 / 6, 2 / 3))
  expect_equal(b$cum.futility, c(1 / 3, 1 / 3))
  expect_equal(unname(b$stage), rbind(c(1 / 6, 1 / 2), c(1 / 3, 0)))
  expect_identical(rownames(b$stage), c("efficacy", "futility"))
  expect_identical(b$names, c("1", "2"))
  expect_equal(pd$xlim, c(1, 2))
  expect_equal(pd$ylim$cumulative, c(0, 18 / 25))
  expect_equal(pd$ylim$stage, c(0, 27 / 50))
})

test_that("simsummary_plot_data places the looks by events or by time", {
  df  <- spd_data()
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96),
                         fut.col = "cox.coef", futility = c(0.2, NA))

  pe <- simsummary_plot_data(res, xaxis = "events")
  expect_identical(pe$xcol, "n.event.mean")
  expect_equal(pe$blocks[[1L]]$x, c(50, 100))
  expect_identical(pe$blocks[[1L]]$names, c("50", "100"))
  expect_equal(pe$xlim, c(50, 100))
  expect_equal(pe$blocks[[1L]]$cum.efficacy, c(1 / 6, 2 / 3))

  pt <- simsummary_plot_data(res, xaxis = "time")
  expect_identical(pt$xcol, "cutoff.mean")
  expect_equal(pt$blocks[[1L]]$x, c(12, 24))
  expect_identical(pt$blocks[[1L]]$names, c("12", "24"))
  expect_equal(pt$xlim, c(12, 24))
})

test_that("simsummary_plot_data separates the blocks", {
  df  <- spd_data()
  dfb <- df
  dfb$logrank.z <- dfb$logrank.z + 0.5
  dp  <- rbind(cbind(df, population = "A"), cbind(dfb, population = "B"))
  resp <- simsummary_fast(dp, eff.col = "logrank.z",
                          efficacy = c(-2.8, -1.96),
                          fut.col = "cox.coef", futility = c(0.2, NA))

  pd <- simsummary_plot_data(resp)
  expect_identical(vapply(pd$blocks, function(b) b$label, ""),
                   c("Population: A", "Population: B"))
  a <- pd$blocks[[1L]]
  b <- pd$blocks[[2L]]
  expect_equal(a$cum.efficacy, c(1 / 6, 2 / 3))
  expect_equal(b$cum.efficacy, c(0, 1 / 3))
  expect_equal(b$cum.futility, c(1 / 3, 1 / 3))
  expect_equal(unname(b$stage), rbind(c(0, 1 / 3), c(1 / 3, 0)))
  expect_equal(pd$ylim$cumulative, c(0, 18 / 25))
  expect_equal(pd$ylim$stage, c(0, 27 / 50))
})

test_that("simsummary_plot_data without futility boundaries", {
  df   <- spd_data()
  resP <- simsummary_fast(df, p.col = "logrank.p", alpha = c(0.0026, 0.025))

  pd <- simsummary_plot_data(resP)
  expect_false(pd$has.futility)
  b <- pd$blocks[[1L]]
  expect_null(b$cum.futility)
  expect_equal(b$cum.efficacy, c(1 / 6, 2 / 3))
  expect_equal(unname(b$stage), rbind(c(1 / 6, 1 / 2)))
  expect_identical(rownames(b$stage), "efficacy")
  expect_equal(pd$ylim$cumulative, c(0, 18 / 25))
  expect_equal(pd$ylim$stage, c(0, 27 / 50))
})

test_that("simsummary_plot_data rejects results that cannot be plotted", {
  df  <- spd_data()
  res <- simsummary_fast(df, eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  expect_error(simsummary_plot_data(res[res$look != "1", ]),
               "keep all their looks")
  expect_error(simsummary_plot_data(res, xaxis = "bars"), "should be one of")
  res_nc <- simsummary_fast(df[, setdiff(names(df), "cutoff")],
                            eff.col = "logrank.z", efficacy = c(-2.8, -1.96))
  expect_error(simsummary_plot_data(res_nc, xaxis = "time"), "cutoff.mean")
  old <- res
  bd <- attr(old, "boundary")
  bd$looks <- NULL
  attr(old, "boundary") <- bd
  expect_error(simsummary_plot_data(old), "earlier version of FastSurvival")
})
