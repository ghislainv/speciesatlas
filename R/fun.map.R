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
    total.taxon.present <- stack()
    taxon.first <- TRUE
    while(taxon.sp[i]==sp){
      if(enough[i]==TRUE){
        nbCategory <- dim(stack(paste0("BIOMOD/",sp.dir[i],"/proj_current/proj_current_",sp.dir[i],".grd")))[[3]]+1##Adding +1 for zeros
        ca <- stack(paste0("BIOMOD/",sp.dir[i],"/proj_current/proj_current_",sp.dir[i],"_ensemble.grd"))[[1]]

        treshold <- seq(0,1000,length.out=nbCategory)[floor(nbCategory/2)+1]
        ca[values(ca)<treshold] <- 0
        ca[values(ca)>=treshold] <- 1

        if(taxon.first==TRUE){ total.taxon.present <- ca
        taxon.first <- FALSE
        } else { total.taxon.present <- total.taxon.present + ca }
      }
      if(i==length(taxon.sp)){sp=0
      } else { i <- i+1 }
    }

    colors <- c("#FFFFFF",colorRampPalette(c("khaki2","orange","red","black"))(max(values(total.taxon.present),na.rm=T)))
    breakpoints <- -0.5:(max(values(total.taxon.present),na.rm=TRUE)+0.5)
    a.arg <- list(at=c(0,max(values(total.taxon.present),na.rm=TRUE)), labels=as.character(c(0,max(values(total.taxon.present),na.rm=TRUE))),cex.axis=1.5)
    l.arg <- list(text="Number of species",side=2, line=0.5, cex=2.5)

    png(paste0("figures/",taxon.names[j],".current.richness.png"),width=650,height=1000)
    plot(total.taxon.present,col=colors,breaks=breakpoints,ext=ext,
         legend.width=1.5,legend.shrink=0.6,legend.mar=7,
         axis.args=a.arg,legend.arg=l.arg,
         axes=FALSE, box=FALSE, zlim=c(0,max(values(total.taxon.present),na.rm=TRUE)))
    dev.off()
    maxvalues <- c(maxvalues,max(values(total.taxon.present),na.rm=TRUE))
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

            # Value used to consider a cell as favorable for the species in the future
            nbCategory <- dim(stack(paste0("BIOMOD/",sp.dir[i],"/proj_current/proj_current_",sp.dir[i],".grd")))[[3]]+1##Adding +1 for zeros
            nbCategoryFuture <- (nbCategory-1)*length(fut.var[[1]])+1

            treshold <- seq(0,1000,length.out=nbCategory)[floor(nbCategory/2)+1]
            tresholdFuture <- seq(0,1000*length(fut.var[[1]]),length.out=nbCategoryFuture)[floor(nbCategoryFuture/2)+1]

            values(caFut)[values(ca)<=treshold] <- 0
            caFut[values(caFut)<tresholdFuture] <- 0
            caFut[values(caFut)>=tresholdFuture] <- 1

            if(taxon.first==TRUE){ total.taxon <- caFut
            taxon.first <- FALSE
            } else { total.taxon <- total.taxon + caFut }
          }
          if(i==length(taxon.sp)){sp=0
          } else { i <- i+1 }
        }

        #
        chang.relatif <- raster(total.taxon)
        diff.absolue <- raster(total.taxon)
        turnover <- raster(total.taxon)

        #Changement relatif au nb esp dans le présent:
        values(chang.relatif) <-((values(total.taxon)-values(total.taxon.present))/values(total.taxon.present))

        #Changement absolue:
        values(diff.absolue) <-abs((values(total.taxon)-values(total.taxon.present)))

        #turn-over
        spGained <- values(total.taxon-total.taxon.present)
        spGained[spGained<0] <- 0

        spLost <- values(total.taxon.present-total.taxon)
        spLost[spLost<0] <- 0

        values(turnover) <- (spGained+spLost)/(values(total.taxon.present)+spGained)

        colors <- c("#FFFFFF",colorRampPalette(c("khaki2","orange","red","black"))(max(values(total.taxon),na.rm=T)))
        breakpoints <- -0.5:(max(values(total.taxon),na.rm=T)+0.5)
        a.arg <- list(at=c(0,max(values(total.taxon),na.rm=T)), labels=as.character(c(0,max(values(total.taxon),na.rm=T))),cex.axis=1.5)
        l.arg <- list(text="Number of species",side=2, line=0.5, cex=2.5)

        png(paste0("figures/",taxon.names[k],".",fut.var[[2]][j],"_",fut.var[[3]][l],".richness.png"),width=650,height=1000)
        plot(total.taxon,col=colors,breaks=breakpoints,ext=ext,
             legend.width=1.5,legend.shrink=0.6,legend.mar=7,
             axis.args=a.arg,legend.arg=l.arg,
             axes=FALSE, box=FALSE, zlim=c(0,max(values(total.taxon),na.rm=T)))
        dev.off()

        colors <- colorRampPalette(c("darkred","yellow","darkgreen"))(100)
        l.arg <- list(text="Perte relative du nombre d'espèces (%)",side=2, line=0.5, cex=1.3)
        png(paste0("figures/",taxon.names[k],".",fut.var[[2]][j],"_",fut.var[[3]][l],".diff.png"),width=650,height=1000)
        plot(chang.relatif*100,col=colors,ext=ext,
             legend.width=1.5,legend.shrink=0.6,legend.mar=7,
             axes=FALSE, box=FALSE,legend.arg=l.arg)
        dev.off()

        colors <- c("#FFFFFF",colorRampPalette(c("khaki2","orange","red","black"))(max(values(diff.absolue),na.rm=T)))
        breakpoints <- c(0.5,1:max(values(diff.absolue),na.rm=T),max(values(diff.absolue),na.rm=T)+0.5)
        a.arg <- list(at=c(0,max(values(diff.absolue),na.rm=T)), labels=as.character(c(0,max(values(diff.absolue),na.rm=T))),cex.axis=1.5)
        l.arg <- list(text="Perte absolue en espèces",side=2, line=0.5, cex=1.8)
        png(paste0("figures/",taxon.names[k],".",fut.var[[2]][j],"_",fut.var[[3]][l],".perte.png"),width=650,height=1000)
        plot(diff.absolue,col=colors,breaks=breakpoints,ext=ext,
             legend.width=1.5,legend.shrink=0.6,legend.mar=7,
             axis.args=a.arg,legend.arg=l.arg,
             axes=FALSE, box=FALSE)
        dev.off()
      }
    }
  }
}
