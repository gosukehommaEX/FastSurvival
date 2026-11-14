# Validation of the switching rules of simdata_fast() (FastSurvival 2.0.0)
# against the multistate model of DesignComplexSurvival.
#
# Run from the package root with the version to be validated installed (for
# example after devtools::install()) or loaded with devtools::load_all():
#   source("inst/validation/validate_switching.R")
# It reads inst/validation/reference/:
#   switching_dcs.csv, switching_dcs_info.txt  analytic values of
#     DesignComplexSurvival (make_switching_reference.R)
#   sim_python.csv, sim_python_two.csv  patient-level simulations in Python
#     of DesignComplexSurvival (see SOURCES.md)
# and writes to inst/validation/results/:
#   validate_switching_summary.md   summary, generated from the results below
#   validate_switching_results.csv  every comparison and self-test
#
# Comparisons, for the scenarios of switching_scenarios.R:
#   (a) Analytic values. The state occupation probabilities and the OS and
#       PFS survival at the latent times, from 10^6 subjects per group
#       without dropout (alternative hypothesis), with
#       z = (p - p0) / sqrt(p0 (1 - p0) / n); and the mean numbers of events
#       and dropouts of OS and PFS per group at the calendar looks, from the
#       simulated trials, with z = (m - mu) / (s / sqrt(R)).
#   (b) Python simulations. The means of the numbers of events and dropouts
#       and of the log-rank Z (positive when the treatment is better) at the
#       calendar looks and at the looks driven by OS deaths, the mean time of
#       those looks, the correlations of Z between looks and between PFS and
#       OS, and the rejection rates of the O'Brien-Fleming design with the
#       Python bounds, with z = (m1 - m2) / sqrt(se1^2 + se2^2).
# Standard deviations are reported but not judged: simdata_fast() fixes the
# group sizes and the numbers enrolled per accrual interval, which the Python
# simulation draws at random. In V12 the Python simulation also has a
# censoring of PFS only (dropout_pfs = 0.01), so its rows that involve PFS are
# INFO.
#
# Verdicts: PASS when |z| <= z_crit, the 1 - 0.01 / (2 m) quantile of the
# standard normal distribution for the m judged comparisons, so that all of
# them are PASS with probability 0.99 when the generators agree; EXPLAINED
# when a reason registered for this script in inst/validation/explained.csv
# applies (none is registered); FAIL otherwise; INFO for rows not judged.
# Self-tests check that the comparisons are not vacuous: data generated with
# one rule changed must give at least one FAIL against the analytic values of
# the scenario: the latent quantities of V7 with switch.clock = "reset", of V9
# with h12.switch.clock = "same", of V10 with switch.aft.factor = 0.6 (the
# factor of DesignComplexSurvival instead of its inverse), and of V11 with
# h01.switch.* left at its default, and the PFS events and dropouts of V12
# with e1.censor.switch = FALSE.

if (!"FastSurvival" %in% loadedNamespaces()) library(FastSurvival)
stopifnot("h1s.hazard" %in% names(formals(FastSurvival::simdata_fast)))
t_start <- Sys.time()

val_dir <- file.path("inst", "validation")
source(file.path(val_dir, "switching_scenarios.R"))
ref_dir <- file.path(val_dir, "reference")
res_dir <- file.path(val_dir, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

explained <- utils::read.csv(file.path(val_dir, "explained.csv"),
                             stringsAsFactors = FALSE)
explained <- explained[explained$script == "validate_switching.R", ]
stopifnot(nrow(explained) == 0L)

read_ref <- function(f, cls) {
  d <- utils::read.csv(file.path(ref_dir, f), stringsAsFactors = FALSE,
                       colClasses = cls)
  for (v in names(cls)) d[[v]][is.na(d[[v]])] <- ""
  d
}
py_cls <- c(look = "character", endpoint = "character", group = "character")
dcs <- read_ref("switching_dcs.csv",
                c(endpoint = "character", group = "character"))
dcs$value <- as.numeric(dcs$value)
py1 <- read_ref("sim_python.csv", py_cls)
py2 <- read_ref("sim_python_two.csv", py_cls)
dcs_info <- readLines(file.path(ref_dir, "switching_dcs_info.txt"))

seed0 <- 20261113L
eps <- c(os = "e2", pfs = "e1")
n_trials <- function(nm) {
  if (nm == "V8") c(R = 4000L, batch = 500L) else c(R = 20000L, batch = 2000L)
}

# ---- Simulated trials: events, dropouts, Z, and look times ----------------
sim_trials <- function(nm, hyp, mod = identity) {
  des <- sw_design(nm)
  rb <- n_trials(nm)
  nb <- rb[["batch"]]
  args <- mod(sw_fs_args(nm, hyp))
  sid <- match(nm, sw_names)
  st <- list()
  add <- function(key, v) st[[key]] <<- c(st[[key]], v)
  for (b in seq_len(rb[["R"]] %/% nb)) {
    d <- do.call(simdata_fast,
                 c(list(nsim = nb, stream = b,
                        seed = seed0 + 10L * sid + (hyp == "null")), args))
    trt <- d$group == 2L
    cutm <- cutoff_fast(d, event.looks = des$events, tte.col = "e2_tte",
                        event.col = "e2_event")
    cutm <- cutm[match(seq_len(nb), as.integer(rownames(cutm))), ,
                 drop = FALSE]
    looks <- list()
    for (k in seq_along(des$cal)) {
      looks[[length(looks) + 1L]] <- list(type = "calendar", k = k,
                                          cut = rep(des$cal[k], nb))
    }
    for (k in seq_along(des$events)) {
      looks[[length(looks) + 1L]] <- list(type = "event", k = k,
                                          cut = unname(cutm[, k]))
      add(paste0("event|", k, "|time"), unname(cutm[, k]))
    }
    for (lk in looks) {
      # As in analysis_fast(), a look not reached uses all data
      cut_i <- lk$cut[d$sim]
      cut_i[is.na(cut_i)] <- Inf
      tag <- paste0(lk$type, "|", lk$k)
      for (ep in names(eps)) {
        tt <- d[[paste0(eps[[ep]], "_tte")]]
        ev <- d[[paste0(eps[[ep]], "_event")]] == 1L
        # On the calendar scale with the sum (accrual + time) of cutoff_fast()
        seen <- is.finite(tt) & d$accrual_time + tt <= cut_i
        m <- rowsum(cbind(ec = seen & ev & !trt, et = seen & ev & trt,
                          dc = seen & !ev & !trt, dt = seen & !ev & trt) * 1,
                    d$sim, reorder = TRUE)
        add(paste0(tag, "|events|", ep, "|control"), m[, "ec"])
        add(paste0(tag, "|events|", ep, "|treatment"), m[, "et"])
        add(paste0(tag, "|dropouts|", ep, "|control"), m[, "dc"])
        add(paste0(tag, "|dropouts|", ep, "|treatment"), m[, "dt"])
        if (lk$type == "event" && ep == "os") {
          reached <- !is.na(lk$cut)
          tot <- (m[, "ec"] + m[, "et"])[reached]
          if (any(tot != des$events[lk$k])) {
            stop(sprintf("%s %s batch %d look %d: the number of OS deaths ",
                         nm, hyp, b, lk$k),
                 "differs from the target in ", sum(tot != des$events[lk$k]),
                 " trials", call. = FALSE)
          }
        }
        de <- data.frame(sim = d$sim, group = d$group,
                         accrual_time = d$accrual_time, tte = tt,
                         event = as.integer(ev))
        cutv <- matrix(lk$cut, ncol = 1L, dimnames = list(seq_len(nb), NULL))
        a <- analysis_fast(de, control = 1, cutoff.looks = cutv)
        a <- a[order(a$sim), ]
        add(paste0(tag, "|z|", ep), -a$logrank.z)
      }
    }
    message(sprintf("  %s %s: batch %d of %d", nm, hyp, b,
                    rb[["R"]] %/% nb))
  }
  st
}

# Summary rows of the simulated trials, with the keys of the Python files
summ_row <- function(nm, hyp, lt, look, time, q, ep, g, x) {
  x <- x[is.finite(x)]
  s <- if (length(x) > 1L) stats::sd(x) else NA_real_
  data.frame(scenario = nm, hypothesis = hyp, look_type = lt, look = look,
             time = time, quantity = q, endpoint = ep, group = g,
             mean = mean(x), sd = s, se = s / sqrt(length(x)), n = length(x),
             stringsAsFactors = FALSE)
}
corr_row <- function(nm, hyp, lt, look, q, ep, z1, z2) {
  ok <- is.finite(z1) & is.finite(z2)
  r <- stats::cor(z1[ok], z2[ok])
  data.frame(scenario = nm, hypothesis = hyp, look_type = lt, look = look,
             time = NA_real_, quantity = q, endpoint = ep, group = "",
             mean = r, sd = NA_real_, se = (1 - r^2) / sqrt(sum(ok)),
             n = sum(ok), stringsAsFactors = FALSE)
}

fs_summary <- function(nm, hyp, st) {
  des <- sw_design(nm)
  rows <- list()
  push <- function(r) rows[[length(rows) + 1L]] <<- r
  for (lt in c("calendar", "event")) {
    nk <- if (lt == "calendar") length(des$cal) else length(des$events)
    for (k in seq_len(nk)) {
      tag <- paste0(lt, "|", k)
      tval <- if (lt == "calendar") des$cal[k] else NA_real_
      if (lt == "event") {
        push(summ_row(nm, hyp, lt, as.character(k), tval, "time", "os", "",
                      st[[paste0(tag, "|time")]]))
      }
      for (ep in names(eps)) {
        for (q in c("events", "dropouts")) {
          for (g in c("control", "treatment")) {
            push(summ_row(nm, hyp, lt, as.character(k), tval, q, ep, g,
                          st[[paste(tag, q, ep, g, sep = "|")]]))
          }
        }
        push(summ_row(nm, hyp, lt, as.character(k), tval, "z", ep, "",
                      st[[paste0(tag, "|z|", ep)]]))
      }
    }
    for (ep in names(eps)) {
      push(corr_row(nm, hyp, lt, "1-2", "z_corr", ep,
                    st[[paste0(lt, "|1|z|", ep)]],
                    st[[paste0(lt, "|2|z|", ep)]]))
    }
  }
  for (j in seq_along(des$cal)) {
    for (k in seq_along(des$cal)) {
      push(corr_row(nm, hyp, "calendar", paste0(j, "-", k), "z_corr_po",
                    "pfs-os", st[[paste0("calendar|", j, "|z|pfs")]],
                    st[[paste0("calendar|", k, "|z|os")]]))
    }
  }
  # Rejections of the O'Brien-Fleming design with the Python bounds
  bnd <- py1[py1$scenario == nm & py1$hypothesis == hyp &
               py1$look_type == "event" & py1$quantity == "bound", ]
  if (nrow(bnd) == 2L) {
    b <- as.numeric(bnd$mean[order(bnd$look)])
    z1 <- st[["event|1|z|os"]]
    z2 <- st[["event|2|z|os"]]
    r1 <- !is.na(z1) & z1 >= b[1L]
    ra <- r1 | (!is.na(z2) & z2 >= b[2L])
    for (q in c("reject_look1", "reject_any")) {
      x <- if (q == "reject_look1") r1 else ra
      p <- mean(x)
      push(data.frame(scenario = nm, hypothesis = hyp, look_type = "event",
                      look = "", time = NA_real_, quantity = q,
                      endpoint = "os", group = "", mean = p, sd = NA_real_,
                      se = sqrt(p * (1 - p) / length(x)), n = length(x),
                      stringsAsFactors = FALSE))
    }
  }
  do.call(rbind, rows)
}

# ---- Latent quantities without dropout -------------------------------------
sim_latent <- function(nm, mod = identity, per_batch = 250000L,
                       n_batch = 4L) {
  des <- sw_design(nm)
  args <- mod(sw_fs_args(nm, "alternative"))
  args$n <- c(per_batch, per_batch)
  args$a.time <- c(0, 1)
  args$a.rate <- NULL
  args$a.prop <- 1
  args$d.hazard <- NULL
  args$d.time <- NULL
  sid <- match(nm, sw_names)
  cnt <- list(control = 0, treatment = 0)
  for (b in seq_len(n_batch)) {
    d <- do.call(simdata_fast, c(list(nsim = 1, stream = b,
                                      seed = seed0 + 10L * sid + 5L), args))
    for (g in names(cnt)) {
      dg <- d[d$group == match(g, names(cnt)), ]
      m <- t(vapply(des$latent, function(t) colSums(sw_states(dg, t)),
                    numeric(7L)))
      cnt[[g]] <- cnt[[g]] + m
    }
  }
  n <- per_batch * n_batch
  rows <- list()
  for (g in names(cnt)) {
    for (q in colnames(cnt[[g]])) {
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = nm, group = g, time = des$latent, quantity = q,
        p = cnt[[g]][, q] / n, n = n, stringsAsFactors = FALSE)
    }
  }
  message(sprintf("  %s: latent quantities", nm))
  do.call(rbind, rows)
}

# ---- Comparisons -----------------------------------------------------------
cmp_latent <- function(lat, ref_nm, source, note = "") {
  r <- dcs[dcs$scenario == ref_nm & dcs$hypothesis == "alternative" &
             grepl("^(p_|surv_)", dcs$quantity), ]
  m <- merge(lat, r, by = c("group", "time", "quantity"),
             suffixes = c("", ".ref"))
  stopifnot(nrow(m) == nrow(lat))
  # A probability of 0 or 1 (up to rounding) has no binomial error: the
  # simulated value must then be equal to it.
  v <- pmin(pmax(m$value, 0), 1)
  se0 <- sqrt(v * (1 - v) / m$n)
  z <- ifelse(se0 > 0, (m$p - m$value) / se0,
              ifelse(abs(m$p - m$value) < 1e-12, 0, Inf))
  data.frame(source = source, scenario = unique(lat$scenario),
             hypothesis = "alternative", look_type = "latent", look = "",
             time = m$time, quantity = m$quantity, endpoint = "",
             group = m$group, fs_mean = m$p, fs_sd = NA_real_, fs_se = se0,
             fs_n = m$n, ref_mean = m$value, ref_sd = NA_real_, ref_se = 0,
             ref_n = NA_real_, z = z, verdict = "", note = note,
             stringsAsFactors = FALSE)
}

cmp_events <- function(fs, ref_nm, source, note = "", endpoints = names(eps)) {
  f <- fs[fs$look_type == "calendar" &
            fs$quantity %in% c("events", "dropouts") &
            fs$endpoint %in% endpoints, ]
  r <- dcs[dcs$scenario == ref_nm & dcs$quantity %in% c("events", "dropouts"), ]
  m <- merge(f, r, by = c("hypothesis", "quantity", "endpoint", "group",
                          "time"), suffixes = c("", ".ref"))
  stopifnot(nrow(m) == nrow(f))
  data.frame(source = source, scenario = m$scenario, hypothesis = m$hypothesis,
             look_type = "calendar", look = m$look, time = m$time,
             quantity = m$quantity, endpoint = m$endpoint, group = m$group,
             fs_mean = m$mean, fs_sd = m$sd, fs_se = m$se, fs_n = m$n,
             ref_mean = m$value, ref_sd = NA_real_, ref_se = 0,
             ref_n = NA_real_, z = (m$mean - m$value) / m$se, verdict = "",
             note = note, stringsAsFactors = FALSE)
}

# Python rows comparable with fs_summary(): sim_python.csv (V6 to V10) and,
# from sim_python_two.csv, the PFS-OS correlations at the calendar looks and,
# for V11 and V12, the looks of the design (calendar months 20 and 30)
py_rows <- function(nm) {
  q1 <- c("events", "dropouts", "z", "z_corr", "time", "reject_look1",
          "reject_any")
  a <- py1[py1$scenario == nm & py1$look_type %in% c("calendar", "event") &
             py1$quantity %in% q1, ]
  b <- py2[py2$scenario == nm & py2$look_type == "calendar" &
             py2$quantity == "z_corr_po", ]
  c2 <- py2[py2$scenario == nm & nm %in% c("V11", "V12") &
              py2$look_type == "design" &
              (py2$quantity %in% c("events", "z") |
                 (py2$quantity == "z_corr" & py2$endpoint == "os-os")), ]
  # c2 has no rows for the scenarios other than V11 and V12
  c2$look_type <- rep("calendar", nrow(c2))
  c2$endpoint[c2$endpoint == "os-os"] <- "os"
  c2$time <- ifelse(c2$quantity == "z_corr", NA_real_,
                    sw_design(nm)$cal[suppressWarnings(as.integer(c2$look))])
  out <- rbind(a, b, c2)
  out[, c("scenario", "hypothesis", "look_type", "look", "quantity",
          "endpoint", "group", "mean", "sd", "se", "n")]
}

cmp_python <- function(fs, nm) {
  py <- py_rows(nm)
  m <- merge(fs, py, by = c("scenario", "hypothesis", "look_type", "look",
                            "quantity", "endpoint", "group"),
             suffixes = c("", ".ref"))
  stopifnot(nrow(m) == nrow(py))
  info <- nm == "V12" & grepl("pfs", m$endpoint)
  data.frame(source = "python", scenario = m$scenario,
             hypothesis = m$hypothesis, look_type = m$look_type,
             look = m$look, time = m$time, quantity = m$quantity,
             endpoint = m$endpoint, group = m$group, fs_mean = m$mean,
             fs_sd = m$sd, fs_se = m$se, fs_n = m$n, ref_mean = m$mean.ref,
             ref_sd = m$sd.ref, ref_se = m$se.ref, ref_n = m$n.ref,
             z = (m$mean - m$mean.ref) / sqrt(m$se^2 + m$se.ref^2),
             verdict = ifelse(info, "INFO", ""),
             note = ifelse(info, "dropout_pfs in the Python simulation", ""),
             stringsAsFactors = FALSE)
}

# ---- Run -------------------------------------------------------------------
res <- list()
fs_all <- list()
lat_all <- list()
for (nm in sw_names) {
  lat_all[[nm]] <- sim_latent(nm)
  res[[length(res) + 1L]] <- cmp_latent(lat_all[[nm]], nm, "dcs_latent")
  fs <- do.call(rbind, lapply(c("alternative", "null"), function(hyp) {
    fs_summary(nm, hyp, sim_trials(nm, hyp))
  }))
  fs_all[[nm]] <- fs
  res[[length(res) + 1L]] <- cmp_events(fs, nm, "dcs_events")
  res[[length(res) + 1L]] <- cmp_python(fs, nm)
}
main <- do.call(rbind, res)

# Self-tests: one rule changed in the arguments of simdata_fast()
set_arg <- function(nm_arg, value) {
  function(a) {
    a[nm_arg] <- list(value)
    a
  }
}
drop_args <- function(a) {
  a$h01.switch.hazard <- NULL
  a$h01.switch.time <- NULL
  a
}
self <- rbind(
  cmp_latent(sim_latent("V7", set_arg("switch.clock", "reset")), "V7",
             "self-test", "V7 generated with switch.clock = \"reset\""),
  cmp_latent(sim_latent("V9", set_arg("h12.switch.clock", "same")), "V9",
             "self-test", "V9 generated with h12.switch.clock = \"same\""),
  cmp_latent(sim_latent("V10", set_arg("switch.aft.factor", 0.6)), "V10",
             "self-test", "V10 generated with switch.aft.factor = 0.6"),
  cmp_latent(sim_latent("V11", drop_args), "V11", "self-test",
             "V11 generated with h01.switch.* at its default"),
  cmp_events(fs_summary("V12", "alternative",
                        sim_trials("V12", "alternative",
                                   set_arg("e1.censor.switch", FALSE))),
             "V12", "self-test", "V12 generated with e1.censor.switch = FALSE",
             endpoints = "pfs"))

judged <- main$verdict != "INFO"
n_judged <- sum(judged)
z_crit <- stats::qnorm(1 - 0.01 / (2 * n_judged))
pass <- function(z) !is.na(z) & abs(z) <= z_crit
main$verdict[judged] <- ifelse(pass(main$z[judged]), "PASS", "FAIL")
self$verdict <- ifelse(pass(self$z), "PASS", "FAIL")
out <- rbind(main, self)
num <- vapply(out, is.numeric, logical(1L))
out[num] <- lapply(out[num], function(x) signif(x, 10))
utils::write.csv(out, file.path(res_dir, "validate_switching_results.csv"),
                 row.names = FALSE)

# ---- Summary ---------------------------------------------------------------
minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
src_lab <- c(dcs_latent = "(a) latent quantities, analytic values",
             dcs_events = "(a) events and dropouts, analytic values",
             python = "(b) Python simulations")
tab <- vapply(names(src_lab), function(s) {
  v <- main$verdict[main$source == s]
  sprintf("| %s | %d | %d | %d | %d | %d |", src_lab[[s]], length(v),
          sum(v == "PASS"), sum(v == "EXPLAINED"), sum(v == "FAIL"),
          sum(v == "INFO"))
}, character(1L))
max_pass <- vapply(names(src_lab), function(s) {
  z <- main$z[main$source == s & main$verdict == "PASS"]
  sprintf("%s %.2f", src_lab[[s]], if (length(z)) max(abs(z)) else NA_real_)
}, character(1L))
st_lines <- vapply(unique(self$note), function(nt) {
  v <- self$verdict[self$note == nt]
  sprintf("- %s: %d of %d comparisons FAIL (largest |z| %.1f).", nt,
          sum(v == "FAIL"), length(v), max(abs(self$z[self$note == nt])))
}, character(1L))
fails <- main[main$verdict == "FAIL", ]
fail_lines <- if (nrow(fails) == 0L) {
  "None."
} else {
  c(paste("| Source | Scenario | Hypothesis | Look | Quantity | Endpoint |",
          "Group | z |"),
    "|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %s | %s | %s %s | %s | %s | %s | %.2f |", fails$source,
            fails$scenario, fails$hypothesis, fails$look_type, fails$look,
            fails$quantity, fails$endpoint, fails$group, fails$z))
}
md <- c(
  "# Validation of the switching rules of simdata_fast()",
  "",
  paste0("Generated by `inst/validation/validate_switching.R` on ",
         format(Sys.Date()), " with ", R.version.string, " and FastSurvival ",
         utils::packageVersion("FastSurvival"), " in ",
         sprintf("%.0f", minutes), " minutes. Do not edit by hand."),
  "",
  paste0("Analytic values: ", dcs_info[1L], " (`reference/switching_dcs.csv`,",
         " computed on ", dcs_info[length(dcs_info)], "). Python simulations:",
         " `reference/sim_python.csv` and `reference/sim_python_two.csv`",
         " (see `reference/SOURCES.md`)."),
  "",
  "| Comparison | Rows | PASS | EXPLAINED | FAIL | INFO |",
  "|---|---|---|---|---|---|",
  tab,
  "",
  paste0("Scenarios: ", paste(sw_names, collapse = ", "),
         " (`switching_scenarios.R`). Simulated trials per hypothesis: 20000",
         " (V8: 4000). Latent quantities: 10^6 subjects per group without",
         " dropout."),
  "",
  sprintf(paste0("Judged comparisons: %d; critical value z_crit = %.3f (the",
                 " 1 - 0.01 / (2 m) quantile). Largest |z| among the PASS",
                 " rows: %s."), n_judged, z_crit,
          paste(max_pass, collapse = "; ")),
  "",
  paste0("INFO: the rows of V12 that involve PFS in the comparison with the",
         " Python simulation, which has a censoring of PFS only",
         " (dropout_pfs = 0.01) that simdata_fast() does not have.",
         " EXPLAINED: no reason is registered for this script."),
  "",
  "FAIL rows:",
  "",
  fail_lines,
  "",
  paste("Self-tests (data generated with one rule changed, against the",
        "analytic values of the scenario):"),
  "",
  st_lines
)
writeLines(md, file.path(res_dir, "validate_switching_summary.md"))

if (any(vapply(unique(self$note), function(nt) {
  !any(self$verdict[self$note == nt] == "FAIL")
}, logical(1L)))) {
  stop("a self-test gave no FAIL: the comparison does not detect the change")
}
