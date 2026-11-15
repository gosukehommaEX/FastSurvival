#' Legend at an explicitly computed keyword position
#'
#' Internal helper: place a legend at a keyword position computed explicitly
#' from par("usr"). Keyword positions passed straight to legend() can misbehave
#' under some multi-panel (mfrow) device states, so the corner coordinates and
#' justification are derived directly and passed as x, y, xjust, yjust.
#'
#' @param pos Legend position: one of the keywords \code{"topright"},
#'   \code{"topleft"}, \code{"bottomright"}, \code{"bottomleft"}, \code{"top"},
#'   \code{"bottom"}, \code{"left"}, \code{"right"}, and \code{"center"}; any
#'   other value is passed to \code{graphics::legend()} unchanged.
#' @param inset Inset from the edges of the plot region, as a fraction of its
#'   width and height.
#' @param ... Further arguments to \code{graphics::legend()}.
#' @return The value of \code{graphics::legend()}.
#' @noRd
legend_fixed <- function(pos, inset = 0.04, ...) {
  usr <- graphics::par("usr")
  dx  <- inset * (usr[2L] - usr[1L])
  dy  <- inset * (usr[4L] - usr[3L])
  xl  <- usr[1L] + dx
  xr  <- usr[2L] - dx
  xc  <- (usr[1L] + usr[2L]) / 2
  yb  <- usr[3L] + dy
  yt  <- usr[4L] - dy
  yc  <- (usr[3L] + usr[4L]) / 2
  spec <- switch(
    pos,
    "topright"    = list(x = xr, y = yt, xj = 1,   yj = 1),
    "topleft"     = list(x = xl, y = yt, xj = 0,   yj = 1),
    "bottomright" = list(x = xr, y = yb, xj = 1,   yj = 0),
    "bottomleft"  = list(x = xl, y = yb, xj = 0,   yj = 0),
    "top"         = list(x = xc, y = yt, xj = 0.5, yj = 1),
    "bottom"      = list(x = xc, y = yb, xj = 0.5, yj = 0),
    "left"        = list(x = xl, y = yc, xj = 0,   yj = 0.5),
    "right"       = list(x = xr, y = yc, xj = 1,   yj = 0.5),
    "center"      = list(x = xc, y = yc, xj = 0.5, yj = 0.5),
    NULL)
  if (is.null(spec)) {
    graphics::legend(pos, ...)
  } else {
    graphics::legend(x = spec$x, y = spec$y, xjust = spec$xj, yjust = spec$yj,
                     ...)
  }
}
