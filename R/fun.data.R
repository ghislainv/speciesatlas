#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurelien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

##=======================
## Data
##=======================

fun.data <- function(df.orig,proj4,run.plots){

  df.sp <- vector("list",length(df.orig))
  sp.names <- NULL
  sp.dir <- NULL
  n.species <- NULL
  taxon.sp <- NULL
  taxon.names <- NULL

  for (sp in 1:length(df.orig)){
    ## Make a SpatialPointsDataFrame object
    coords <- cbind(df.orig[[sp]]$Long,df.orig[[sp]]$Lat)
    df.sp[[sp]] <- SpatialPointsDataFrame(coords,data=as.data.frame(df.orig[[sp]]),proj4string=CRS("+init=epsg:4326"))
    ## Reproject as raster projection
    df.sp[[sp]] <- spTransform(df.sp[[sp]],proj4)
    ## Species
    df.sp[[sp]]$Species <- gsub("\\s$","",df.sp[[sp]]$Species) # Remove final space in names
    sp.names <- c(sp.names,levels(as.factor(df.sp[[sp]]$Species))) # Sorted in alphabetical order
    sp.dir <- c(sp.dir,gsub(" ",".",levels(as.factor(df.sp[[sp]]$Species))))
    n.species <- c(n.species,length(levels(as.factor(df.sp[[sp]]$Species))))
    taxon.sp <- c(taxon.sp,rep(sp,n.species[sp]))
    taxon.names <- c(taxon.names,df.orig[[sp]]$Taxo.group)
    ## Figures directory
    dir.create(paste0("figures/",taxon.names[sp]),recursive=TRUE,showWarnings=FALSE)
    ## initializing empty rda for every taxon only if going to run plots
    if(run.plots){
      SDA.whole.fut <- data.frame(species=character(),area.pres=numeric(),area.fut=numeric(),perc.chang=numeric())
      save(list=c("SDA.whole.fut"),file=paste0("figures/",taxon.names[sp],"/plotting.rda"))
    }
  }

  dir.create("BIOMOD",recursive=TRUE,showWarnings=FALSE)
  setwd("BIOMOD")
  return(list(df.sp,sp.names,sp.dir,taxon.sp,taxon.names))
}
