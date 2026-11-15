#' Draw a single scenario panel
#'
#' @param label Title of the panel.
#' @param ev The output of \code{scenario_eval()}.
#' @param hr_max Upper limit of the hazard-ratio axis.
#' @param col_control,col_treat,col_hr Colours of the control and treatment
#'   survival curves and of the hazard ratio.
#' @param lty_control,lty_treat,lty_hr Their line types.
#' @param lwd_control,lwd_treat,lwd_hr Their line widths.
#' @param base_size Scale of the text.
#' @param xlab,ylab Axis labels.
#' @param legend_pos Legend position passed to \code{legend_fixed()}.
#' @param add_legend Whether to draw the legend.
#' @return Called to draw; the value is not used.
#' @noRd
draw_scenario_panel <- function(label, ev, hr_max,
                                col_control, col_treat, col_hr,
                                lty_control, lty_treat, lty_hr,
                                lwd_control, lwd_treat, lwd_hr,
                                base_size, xlab, ylab, legend_pos,
                                add_legend) {
  graphics::plot.default(ev$times, ev$s_c, type = "l", col = col_control,
                         lwd = lwd_control, lty = lty_control,
                         xaxs = "i", yaxs = "i",
                         ylim = c(0, 1), xlab = xlab,
                         ylab = ylab,
                         main = label, cex.main = 0.95 * base_size)
  if (isTRUE(ev$two)) {
    graphics::lines(ev$times, ev$s_t, col = col_treat,
                    lwd = lwd_treat, lty = lty_treat)
    graphics::par(new = TRUE)
    graphics::plot.default(ev$times, ev$hr / hr_max, type = "l",
                           col = col_hr, lwd = lwd_hr, lty = lty_hr,
                           xaxs = "i", yaxs = "i",
                           ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    graphics::abline(h = 1 / hr_max, col = "grey60", lty = 3)
    hr_ticks <- pretty(c(0, hr_max), n = 4)
    hr_ticks <- hr_ticks[hr_ticks <= hr_max]
    graphics::axis(4, at = hr_ticks / hr_max, labels = hr_ticks)
    graphics::mtext("Hazard ratio", side = 4, line = 2.5, cex = base_size)
  }
  if (add_legend) {
    if (isTRUE(ev$two)) {
      leg  <- c("Control", "Experimental", "Hazard ratio")
      cols <- c(col_control, col_treat, col_hr)
      lwds <- c(lwd_control, lwd_treat, lwd_hr)
      ltys <- c(lty_control, lty_treat, lty_hr)
    } else {
      leg  <- "Survival"
      cols <- col_control
      lwds <- lwd_control
      ltys <- lty_control
    }
    legend_fixed(legend_pos, bty = "n", cex = 0.75 * base_size,
                 seg.len = 2.5, legend = leg, col = cols,
                 lwd = lwds, lty = ltys)
  }
}
