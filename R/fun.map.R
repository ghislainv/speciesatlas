#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurelien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

# ==================
# Creating Maps
# ==================

fun.map <- function(sp.dir,ext,enough,taxon.names,taxon.sp,fut.var){

  # CURRENT
  # with the distribution maps obtained
  i=1
  first <- TRUE
  col <- c("#CD3333","#7FFF00","#6495ED","#FF7F50","#00CDCD","#9A32CD")
  maxvalues <- NULL
  for (j in 1:length(unique(taxon.sp))) {
    sp=taxon.sp[i]
    total.taxon <- stack()
    taxon.first <- TRUE
    while(taxon.sp[i]==sp){
      if(enough[i]==TRUE){
        ca <- stack(paste0("BIOMOD/",sp.dir[i],"/proj_current/proj_current_",sp.dir[i],"_ensemble.grd"))[[1]]
        ca[values(ca)<600] <- 0
        ca[values(ca)>=600] <- 1
        if(taxon.first==TRUE){ total.taxon <- ca
        taxon.first <- FALSE
        } else { total.taxon <- total.taxon + ca }
      }
      if(i==length(taxon.sp)){sp=0
      } else { i <- i+1 }
    }

    gcolors <- colorRampPalette(c("#F2F2F2",col[j]))
    colors <- gcolors((max(values(total.taxon),na.rm=TRUE)+1))
    breakpoints <- -0.5:(max(values(total.taxon),na.rm=TRUE)+0.5)
    a.arg <- list(at=c(0,max(values(total.taxon),na.rm=TRUE)), labels=as.character(c(0,max(values(total.taxon),na.rm=TRUE))),cex.axis=1.5)
    l.arg <- list(text="Number of species",side=2, line=0.5, cex=2.5)

    png(paste0("figures/",taxon.names[j],".current.richness.png"),width=650,height=1000)
    plot(total.taxon,col=colors,breaks=breakpoints,ext=ext,
         legend.width=1.5,legend.shrink=0.6,legend.mar=7,
         axis.args=a.arg,legend.arg=l.arg,
         axes=FALSE, box=FALSE, zlim=c(0,max(values(total.taxon),na.rm=TRUE)))
    dev.off()
    maxvalues <- c(maxvalues,max(values(total.taxon),na.rm=TRUE))
  }



  # FUTURE
  for (j in 1:length(fut.var[[2]])) {
    for (l in 1:length(fut.var[[3]])) {
      i=1
      first <- TRUE
      for (k in 1:length(unique(taxon.sp))) {
        sp=taxon.sp[i]
        total.taxon <- stack()
        taxon.first <- TRUE
        while(taxon.sp[i]==sp){
          if(enough[i]==TRUE){
            caS <- stack()
            for (mc in 1:length(fut.var[[1]])) {
              name.f <- paste0("/proj_",fut.var[[1]][mc],"_",fut.var[[2]][j],"_",fut.var[[3]][l])
              pred.f <- stack(paste0("BIOMOD/",sp.dir[i],name.f,name.f,"_",sp.dir[i],"_ensemble.grd"))
              caS <- addLayer(caS,pred.f[[1]])
            }
            # Sum and binarize
            ca <- stack(paste0("BIOMOD/",sp.dir[i],"/proj_current/proj_current_",sp.dir[i],"_ensemble.grd"))[[1]]
            caFut <- stack()
            caFut <- sum(caS)
            values(caFut)[values(ca)<=600] <- 0
            caFut[values(caFut)<1500] <- 0
            caFut[values(caFut)>=1500] <- 1
            if(taxon.first==TRUE){ total.taxon <- caFut
            taxon.first <- FALSE
            } else { total.taxon <- total.taxon + caFut }
          }
          if(i==length(taxon.sp)){sp=0
          } else { i <- i+1 }
        }

        gcolors <- colorRampPalette(c("#F2F2F2",col[k]))
        colors <- gcolors(maxvalues[k]+1)
        breakpoints <- -0.5:(maxvalues[k]+0.5)
        a.arg <- list(at=c(0,maxvalues[k]), labels=as.character(c(0,maxvalues[k])),cex.axis=1.5)
        l.arg <- list(text="Number of species",side=2, line=0.5, cex=2.5)

        png(paste0("figures/",taxon.names[k],".",fut.var[[2]][j],"_",fut.var[[3]][l],".richness.png"),width=650,height=1000)
        plot(total.taxon,col=colors,breaks=breakpoints,ext=ext,
             legend.width=1.5,legend.shrink=0.6,legend.mar=7,
             axis.args=a.arg,legend.arg=l.arg,
             axes=FALSE, box=FALSE, zlim=c(0,maxvalues[k]))
        dev.off()
      }
    }
  }
}
