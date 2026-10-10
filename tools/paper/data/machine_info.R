# Helper sourced at the top of every script in tools/paper/data/.
#
# It stops unless the installed FastSurvival is the version described in the
# article, and writes the computing environment (CPU, cores, memory, operating
# system, R and package versions) to machine_info.csv in the output folder for
# the computational details of the article. Run the scripts from the package
# root (FastSurvival.Rproj), where they write to tools/paper/output, or from
# the article folder of the supplementary material, where they write to data.
# With options(paper.smoke = TRUE) the scripts run with few simulated trials
# and write to a temporary folder, as a quick check before a full run.

paper_fs_version <- "1.2.0"

if (packageVersion("FastSurvival") != paper_fs_version) {
  stop("The paper scripts expect FastSurvival ", paper_fs_version,
       ", but version ", packageVersion("FastSurvival"), " is installed. ",
       "Install it from CRAN with install.packages(\"FastSurvival\", ",
       "type = \"source\") and restart R.", call. = FALSE)
}

paper_smoke <- isTRUE(getOption("paper.smoke"))
paper_out_dir <- if (paper_smoke) {
  file.path(tempdir(), "paper_smoke")
} else if (dir.exists(file.path("tools", "paper", "data"))) {
  file.path("tools", "paper", "output")
} else {
  "data"
}
dir.create(paper_out_dir, showWarnings = FALSE, recursive = TRUE)
if (paper_smoke) message("Smoke test: results are written to ", paper_out_dir)

# Number of simulated trials (or of timed runs): full in a full run, smoke in a
# smoke test.
paper_n <- function(full, smoke) if (paper_smoke) smoke else full

local({
  # First non-empty line of a system command, or NA if it fails.
  run_cmd <- function(cmd, args) {
    out <- tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE,
                                             stderr = FALSE)),
                    error = function(e) character(0))
    out <- trimws(out)
    out <- out[nzchar(out)]
    if (length(out) == 0L) NA_character_ else out[1L]
  }
  sys <- Sys.info()[["sysname"]]
  if (sys == "Windows") {
    ps <- function(expr) {
      run_cmd("powershell", c("-NoProfile", "-Command", shQuote(expr)))
    }
    cpu <- ps("(Get-CimInstance Win32_Processor | Select-Object -First 1).Name")
    mem <- suppressWarnings(as.numeric(
      ps("(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory")))
  } else if (sys == "Darwin") {
    cpu <- run_cmd("sysctl", c("-n", "machdep.cpu.brand_string"))
    mem <- suppressWarnings(as.numeric(
      run_cmd("sysctl", c("-n", "hw.memsize"))))
  } else {
    cpuinfo <- tryCatch(readLines("/proc/cpuinfo"),
                        error = function(e) character(0))
    meminfo <- tryCatch(readLines("/proc/meminfo"),
                        error = function(e) character(0))
    cpu <- sub(".*:[[:space:]]*", "", grep("^model name", cpuinfo,
                                           value = TRUE)[1L])
    mem <- 1024 * suppressWarnings(as.numeric(gsub(
      "[^0-9]", "", grep("^MemTotal", meminfo, value = TRUE)[1L])))
  }

  pkgs <- c("FastSurvival", "simtrial", "TrialSimulator", "gsDesign",
            "survival", "survRM2", "nph", "nphRCT", "survAH",
            "microbenchmark", "mvtnorm", "dqrng", "Rcpp", "future")
  # Installed versions, read from the DESCRIPTION files without loading.
  ver <- vapply(pkgs, function(p) {
    v <- suppressWarnings(utils::packageDescription(p, fields = "Version"))
    if (is.na(v)) "not installed" else v
  }, character(1))

  info <- data.frame(
    item  = c("run_at", "time_zone", "r_version", "os", "cpu",
              "logical_cores", "physical_cores", "memory_gb",
              paste0("pkg_", pkgs)),
    value = c(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Sys.timezone(),
              R.version.string, utils::osVersion, cpu,
              parallel::detectCores(logical = TRUE),
              parallel::detectCores(logical = FALSE),
              round(mem / 2^30, 1), ver),
    stringsAsFactors = FALSE
  )
  utils::write.csv(info, file.path(paper_out_dir, "machine_info.csv"),
                   row.names = FALSE)
})
