#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurélien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

# ==================
# Automatic zoom
# ==================

fun.extent <- function(p,s) {
  width <- xmax(p)-xmin(p)
  height <- ymax(p)-ymin(p)
  if (width<200000 & height<500000) {
    zoom <- TRUE
    x1 <- floor(xmin(p)-50000)
    x2 <- ceiling(xmax(p)+50000)
    y1 <- floor(ymin(p)-100000)
    y2 <- ceiling(ymax(p)+100000)
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