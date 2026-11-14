# Analytic values of DesignComplexSurvival for the validation of the switching
# rules of simdata_fast() (validate_switching.R): for the scenarios of
# switching_scenarios.R, the state occupation probabilities and the OS and PFS
# survival without dropout (ms_prob() and ms_surv(), alternative hypothesis,
# both groups), and the expected numbers of events and dropouts of OS and PFS
# at the calendar looks (ms_events(), both hypotheses). V12 is computed
# without the censoring of PFS only (dropout_pfs), which simdata_fast() does
# not have.
#
# DesignComplexSurvival (0.4.0 or later) is not a dependency of FastSurvival
# and is needed only here. Run from the package root of FastSurvival, either
# with DesignComplexSurvival installed or with the path of its source folder
# in 'dcs_path' (loaded with pkgload::load_all()), for example:
#   dcs_path <- "../../04_DesignComplexSurvival/Rpackage/DesignComplexSurvival"
#   source("inst/validation/make_switching_reference.R")
# It writes inst/validation/reference/switching_dcs.csv and
# switching_dcs_info.txt; validate_switching.R reads them.

if (!exists("dcs_path")) dcs_path <- NULL
if (is.null(dcs_path)) {
  library(DesignComplexSurvival)
  dcs_source <- "installed package"
} else {
  pkgload::load_all(dcs_path, export_all = FALSE, helpers = FALSE,
                    attach_testthat = FALSE, quiet = TRUE)
  # Only the name of the folder is recorded, not its location on this computer
  dcs_source <- paste("source folder", basename(normalizePath(dcs_path)))
}
dcs_version <- as.character(utils::packageVersion("DesignComplexSurvival"))
stopifnot(utils::compareVersion(dcs_version, "0.4.0") >= 0)

val_dir <- file.path("inst", "validation")
source(file.path(val_dir, "switching_scenarios.R"))
ref_dir <- file.path(val_dir, "reference")
dir.create(ref_dir, showWarnings = FALSE)

rows <- list()
add_rows <- function(scenario, hyp, quantity, endpoint, group, time, value) {
  rows[[length(rows) + 1L]] <<- data.frame(
    scenario = scenario, hypothesis = hyp, quantity = quantity,
    endpoint = endpoint, group = group, time = time, value = value,
    stringsAsFactors = FALSE)
}

for (nm in sw_names) {
  des <- sw_design(nm)
  for (hyp in c("alternative", "null")) {
    m <- sw_dcs_scenario(nm, hyp)
    if (hyp == "alternative") {
      for (g in c("control", "treatment")) {
        tt <- des$latent
        p <- ms_prob(m[[g]], tt)
        pB <- if ("progression_free_switched" %in% names(p)) {
          p$progression_free_switched
        } else {
          rep(0, length(tt))
        }
        st <- list(p_A = p$progression_free, p_B = pB, p_C = p$progressed,
                   p_D = p$switched, p_E = p$dead,
                   surv_os = ms_surv(m[[g]], tt, "os")$surv,
                   surv_pfs = ms_surv(m[[g]], tt, "pfs")$surv)
        for (q in names(st)) add_rows(nm, hyp, q, "", g, tt, st[[q]])
      }
    }
    for (ep in c("os", "pfs")) {
      ev <- ms_events(m$scenario, des$n, des$cal, ep)
      for (g in c("control", "treatment")) {
        add_rows(nm, hyp, "events", ep, g, des$cal,
                 ev[[paste0("events_", g)]])
        add_rows(nm, hyp, "dropouts", ep, g, des$cal,
                 ev[[paste0("dropouts_", g)]])
      }
    }
  }
}

out <- do.call(rbind, rows)
out$value <- sprintf("%.17g", out$value)
utils::write.csv(out, file.path(ref_dir, "switching_dcs.csv"),
                 row.names = FALSE)
writeLines(c(paste("DesignComplexSurvival", dcs_version, "from the",
                   dcs_source),
             R.version.string,
             paste("rows", nrow(out)),
             format(Sys.time(), tz = "UTC", usetz = TRUE)),
           file.path(ref_dir, "switching_dcs_info.txt"))
