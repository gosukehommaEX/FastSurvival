#' Print Method for Sequential Analysis Summaries
#'
#' @description
#' Formats an object returned by \code{\link{simsummary_fast}} for display. The
#' boundary settings are shown first, then a per-population table whose look rows
#' carry the stage-wise stopping probabilities and analysis timing and whose
#' \code{overall} row carries the totals.
#'
#' @details
#' The probability columns are printed as proportions rounded to \code{digits}
#' significant digits, and the timing columns (mean enrolled, mean events, mean
#' dropouts, mean pipeline, mean calendar time) are rounded for readability. When
#' the object holds more than one
#' population, each population is printed as its own block with a heading. The
#' underlying object is an ordinary data frame, so the unrounded values remain
#' available by subsetting it directly.
#'
#' @param x An object of class \code{"simsummary_fast"} from
#'   \code{\link{simsummary_fast}}.
#' @param digits A single positive integer, the number of significant digits used
#'   when rounding the printed values. Defaults to 4.
#' @param ... Further arguments, currently ignored.
#'
#' @return The object \code{x}, invisibly.
#'
#' @seealso \code{\link{simsummary_fast}}.
#'
#' @export
print.simsummary_fast <- function(x, digits = 4, ...) {
  bd   <- attr(x, "boundary")
  nsim <- attr(x, "nsim")

  cat("Sequential analysis summary (simsummary_fast)\n")
  if (!is.null(nsim)) {
    if (length(nsim) == 1L) {
      cat("  Simulations: ", nsim, "\n", sep = "")
    } else {
      cat("  Simulations: ", paste(nsim, collapse = ", "),
          " (per population)\n", sep = "")
    }
  }

  # Boundary description
  if (!is.null(bd)) {
    if (identical(bd$mode, "z")) {
      cat("  Mode: Z-scale boundaries, direction = ", bd$direction, "\n",
          sep = "")
      cat("    efficacy [", bd$eff.col, "]: ",
          paste(format(bd$efficacy, digits = digits), collapse = ", "),
          "\n", sep = "")
      if (!is.null(bd$futility)) {
        cat("    futility [", bd$fut.col, "]: ",
            paste(format(bd$futility, digits = digits), collapse = ", "),
            "\n", sep = "")
      }
    } else {
      cat("  Mode: p-value boundaries\n")
      cat("    alpha [", bd$p.col, "]: ",
          paste(format(bd$alpha, digits = digits), collapse = ", "),
          "\n", sep = "")
    }
  }
  cat("\n")

  # Round numeric columns for display, keep the data frame structure
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  num_cols <- vapply(df, is.numeric, logical(1L))
  df[num_cols] <- lapply(df[num_cols], function(v) round(v, digits))

  pops <- unique(df$population)
  single_pop <- length(pops) == 1L

  for (p in pops) {
    block <- df[df$population == p, , drop = FALSE]
    # Drop the constant population column from the printed block
    block$population <- NULL
    rownames(block) <- NULL
    if (!single_pop) {
      cat("Population: ", p, "\n", sep = "")
    }
    print(block, row.names = FALSE)
    cat("\n")
  }

  invisible(x)
}
