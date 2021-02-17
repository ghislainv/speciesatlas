#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurelien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

# ==================
# Automatic zoom
# ==================

fun.extent <- function(p,s) {
  width <- xmax(p)-xmin(p)
  height <- ymax(p)-ymin(p)
  xlim <- (xmax(s)-xmin(s))/4
  ylim <- (ymax(s)-ymin(s))/3

  if (width<xlim & height<ylim) {
    zoom <- TRUE
    x1 <- floor(xmin(p)-xlim/4)
    x2 <- ceiling(xmax(p)+xlim/4)
    y1 <- floor(ymin(p)-ylim/5)
    y2 <- ceiling(ymax(p)+ylim/5)
    e.map <- c(max(x1,xmin(s)),min(x2,xmax(s)),
               max(y1,ymin(s)),min(y2,ymax(s)))
    r.mar <- 4
  }  else {
    zoom <- FALSE
    e.map <- extent(s)
    r.mar <- 1
  }
  return(list(zoom=zoom,e.map=e.map,r.mar=r.mar))
}
