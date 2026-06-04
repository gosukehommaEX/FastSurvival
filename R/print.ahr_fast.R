#' Print an ahr_fast object
#'
#' @param x an object of class \code{"ahr_fast"}
#' @param digits number of significant digits to print
#' @param ... further arguments (currently ignored)
#' @return \code{x}, invisibly
#' @method print ahr_fast
#' @export
print.ahr_fast <- function(x, digits = 4, ...) {
  cat("\nKalbfleisch-Prentice average hazard ratio\n")
  cat("Comparison group '", as.character(x$groups[2]),
      "' vs reference '", as.character(x$groups[1]),
      "' over [0, ", format(x$tau, digits = digits), "]\n\n", sep = "")

  cl <- round(100 * x$conf.level, 1)
  cat("  AHR          : ", format(x$ahr, digits = digits), "\n", sep = "")
  cat("  ", cl, "% CI      : ", format(x$lower, digits = digits), ", ",
      format(x$upper, digits = digits), "\n", sep = "")
  cat("  log(AHR)     : ", format(x$log.ahr, digits = digits),
      "  (se ", format(x$se.loghr, digits = digits), ")\n", sep = "")
  cat("  Z (H0: AHR = ", format(x$null.ahr, digits = digits), ") : ",
      format(x$z, digits = digits), "  [theta scale]\n", sep = "")
  cat("  p-value      : ", format.pval(x$p.value, digits = digits),
      "  [theta scale]\n", sep = "")
  cat("  (log scale: Z = ", format(x$z.loghr, digits = digits), ", p = ",
      format.pval(x$p.value.loghr, digits = digits), ")\n", sep = "")

  cat("\nGroup shares of total hazard (theta):\n")
  tab <- rbind(theta = round(x$theta, digits), n = x$n)
  print(tab)

  cat("\n")
  invisible(x)
}
