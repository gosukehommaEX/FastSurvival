# Scenarios of the validation of the switching rules of simdata_fast(),
# sourced by make_switching_reference.R and validate_switching.R.
#
# The names and parameters are those of the patient-level simulations of the
# DesignComplexSurvival package (version 0.4.0, inst/validation/python/
# simulate_ms.py and simulate_two.py): scenario C of its Phase 3 plan with
# switching after progression (V6, reset clock), the forward clock (V7), the
# post-switch death hazard measured from the switch (V9, forward clock), an
# accelerated failure time effect of the switch (V10, acceleration factor 0.6
# of DesignComplexSurvival, that is, a factor 1 / 0.6 on the remaining time),
# switching before progression (V11), and the censoring of PFS at that switch
# (V12); and Example 1 of Luo et al. (2019), treatment discontinuation with the
# forward clock (V8). In V12 the Python simulation also has a censoring of PFS
# only (dropout_pfs = 0.01), which simdata_fast() does not have. Under the
# null hypothesis both groups follow the control model.
#
# Luo, X., Mao, X., Chen, X., Qiu, J., Bai, S., & Quan, H. (2019). Design and
# monitoring of survival trials in complex scenarios. Statistics in Medicine,
# 38(2), 192-209.

sw_names <- c("V6", "V7", "V8", "V9", "V10", "V11", "V12")

# Total sample size, calendar looks, targets of the looks driven by OS deaths,
# and the times (since randomization) of the latent quantities
sw_design <- function(name) {
  if (name == "V8") {
    list(n = 4600, cal = c(30, 40), events = c(361, 722),
         latent = c(6, 12, 24, 36))
  } else {
    list(n = 500, cal = c(20, 30), events = c(150, 250),
         latent = c(2, 5, 10, 20, 30))
  }
}

# Parameters of V8 (per month)
sw_v8 <- function() {
  list(tm = c(0, 6, 12, 24),
       lam0 = c(0.09, 0.09, 0.07, 0.05) / 12,
       drop = c(0, 0, 0.015, 0.015) / 12,
       disc = c(0.04, 0.04, 0.03, 0.03) / 12,
       hr = c(1, 0.75, 0.75, 0.75),
       rate = c(100, 200, 300, 300, rep(400, 8), 300, 200))
}

# Arguments of simdata_fast() without nsim, seed, and stream
sw_fs_args <- function(name, hyp) {
  stopifnot(name %in% sw_names, hyp %in% c("alternative", "null"))
  if (name == "V8") {
    p <- sw_v8()
    tm <- c(p$tm, Inf)
    ctl <- list(h01 = 0, h02 = p$lam0)
    trt <- if (hyp == "null") ctl else list(h01 = p$disc, h02 = p$hr * p$lam0)
    return(list(n = c(2300, 2300), a.time = 0:14, a.rate = p$rate,
                h01.hazard = list(ctl$h01, trt$h01), h01.time = tm,
                h02.hazard = list(ctl$h02, trt$h02), h02.time = tm,
                h12.hazard = p$lam0, h12.time = tm,
                d.hazard = p$drop, d.time = tm, switch.clock = "forward"))
  }
  pre <- name %in% c("V11", "V12")
  ctl <- list(h01 = c(0.10, 0.08), h02 = 0.01, h12 = c(0.05, 0.06),
              sp = 0.3, h1s = 0.04, h0s = if (pre) 0.02 else 0,
              h01s = if (pre) c(0.07, 0.06), h02s = if (pre) 0.008)
  trt <- list(h01 = c(0.06, 0.05), h02 = 0.008, h12 = c(0.045, 0.05),
              sp = 0, h1s = 0, h0s = 0, h01s = NULL, h02s = NULL)
  if (hyp == "null") trt <- ctl
  args <- list(n = c(250, 250), a.time = c(0, 6, 18), a.prop = c(0.3, 0.7),
               h01.hazard = list(ctl$h01, trt$h01), h01.time = c(0, 4, Inf),
               h02.hazard = list(ctl$h02, trt$h02),
               h12.hazard = list(ctl$h12, trt$h12), h12.time = c(0, 3, Inf),
               switch.prop = list(ctl$sp, trt$sp),
               h1s.hazard = list(ctl$h1s, trt$h1s), d.hazard = 0.005,
               switch.clock = if (name %in% c("V7", "V9")) "forward" else
                 "reset",
               h12.switch.clock = if (name == "V9") "switch" else "same")
  if (name == "V10") {
    args$switch.aft.factor <- 1 / 0.6
  } else {
    args$h12.switch.hazard <- c(0.03, 0.035)
    args$h12.switch.time <- c(0, 3, Inf)
  }
  if (pre) {
    args$h0s.hazard <- list(ctl$h0s, trt$h0s)
    args$h01.switch.hazard <- list(ctl$h01s, trt$h01s)
    args$h01.switch.time <- c(0, 4, Inf)
    args$h02.switch.hazard <- list(ctl$h02s, trt$h02s)
    args$e1.censor.switch <- name == "V12"
  }
  args
}

# Models and scenario of DesignComplexSurvival (its functions ms_model(),
# ms_accrual(), and ms_scenario() must be available); V12 without dropout_pfs
sw_dcs_scenario <- function(name, hyp) {
  stopifnot(name %in% sw_names, hyp %in% c("alternative", "null"))
  if (name == "V8") {
    p <- sw_v8()
    ctl <- ms_model(time = p$tm, death = p$lam0, dropout = p$drop)
    trt <- ms_model(time = p$tm, prog = p$disc, death = p$hr * p$lam0,
                    death_prog = p$lam0, time_post = p$tm, dropout = p$drop,
                    clock = "forward")
    acc <- ms_accrual(time = 0:13, rate = p$rate)
    pfs_switch <- "event"
  } else {
    pre <- name %in% c("V11", "V12")
    clock <- if (name %in% c("V7", "V9")) "forward" else "reset"
    ctl <- ms_model(time = c(0, 4), prog = c(0.10, 0.08), death = 0.01,
                    death_prog = c(0.05, 0.06), switch_prob = 0.3,
                    switch_rate = 0.04,
                    death_switch = if (name != "V10") c(0.03, 0.035),
                    dropout = 0.005, time_post = c(0, 3), clock = clock,
                    switch_clock = if (name == "V9") "switch" else "same",
                    switch_af = if (name == "V10") 0.6,
                    switch_pre = if (pre) 0.02 else 0,
                    prog_switch = if (pre) c(0.07, 0.06),
                    death_switch_pre = if (pre) 0.008)
    trt <- ms_model(time = c(0, 4), prog = c(0.06, 0.05), death = 0.008,
                    death_prog = c(0.045, 0.05), dropout = 0.005,
                    time_post = c(0, 3), clock = clock)
    acc <- ms_accrual(time = c(0, 6), prop = c(0.3, 0.7), duration = 18)
    pfs_switch <- if (name == "V12") "censor" else "event"
  }
  if (hyp == "null") trt <- ctl
  list(control = ctl, treatment = trt,
       scenario = ms_scenario(ctl, trt, allocation = 1, accrual = acc,
                              pfs_switch = pfs_switch))
}

# State occupation of the data at time t (latent columns): A event-free
# without switch, B event-free after a switch before the intermediate event,
# C after the intermediate event without switch, D switched after the
# intermediate event (or after an intermediate event that followed a switch),
# E after the terminal event; and the OS and PFS survival.
sw_states <- function(d, t) {
  sw_t <- d$switched == 1L & d$switch_time <= t
  sw_t[is.na(sw_t)] <- FALSE
  free <- d$e1_surv_time > t
  dead <- d$e2_surv_time <= t
  cbind(p_A = free & !sw_t, p_B = free & sw_t, p_C = !free & !dead & !sw_t,
        p_D = !free & !dead & sw_t, p_E = dead, surv_os = !dead,
        surv_pfs = free)
}
