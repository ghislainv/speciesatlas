#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurelien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

# ==================
# Main function
# ==================

fun.species <- function(i,run.models,run.plots,run.taxo,model.var,environ,future,fut.var,maxent.path,out.type,taxon.sp,taxon.names,sp.dir,sp.names,df.sp,s){

  ## Simplifications
  sp <- taxon.sp[i]
  name <- taxon.names[sp]
  spdir <- sp.dir[i]
  spname <- sp.names[i]
  path <- paste0("../figures/",name,"/",spdir)

  ##===============================
  ## Select data for target species
  df.tsp <- df.sp[[sp]][df.sp[[sp]]$Species==spname,]

  ##===================
  ## Remove duplicates
  cell.pres <- cellFromXY(s,df.tsp)

  if(!all(is.na(cell.pres))){
    ## Spatial points at center of each raster cell with presences
    list.cell.pres <- sort(unique(cell.pres))
    cell.pres.sp <- xyFromCell(s,list.cell.pres,spatial=TRUE)

    ## Build data-set for presence
    enough <- TRUE
    dir.create(path,recursive=TRUE,showWarnings=FALSE)
    d.presence <- as.data.frame(extract(s,cell.pres.sp))
    Coords.presence <- coordinates(cell.pres.sp)
    colnames(Coords.presence) <- c("x","y")
    data.xy <- cbind(d.presence,Coords.presence)

    ## Uncomplete data points
    wcomp <- which(complete.cases(data.xy))
  } else {wcomp <- numeric(0)}


  ## Check if there's enough observations
  if (length(wcomp)<=10){
    enough <- FALSE
    if (length(wcomp)>0){
      p <- SpatialPoints(coords = Coords.presence)
    }
    zoom <-FALSE
    npix<- length(wcomp)

  }else{
    ## Transform as a SpatialPointsDataFrame and SpatialPoints (for presence only)
    d <- SpatialPointsDataFrame(coords=Coords.presence[wcomp,], data=data.xy[wcomp,],
                                proj4string=crs(s))
    p <- SpatialPoints(d, proj4string=crs(d)) ## This is used for presence-only data



    ## Extent for species distribution area maps
    ext <- fun.extent(p,s)
    zoom <- ext$zoom
    e.map <- ext$e.map
    r.mar <- ext$r.mar

    ## Number of 1km pixels with at least one presence
    npix <- nrow(d)

    ##==========
    ## Modelling

    dir.create(spdir,recursive=TRUE,showWarnings=FALSE)

    if(run.models){
      Biomod=fun.models.run(name,spdir,p,s,spname,model.var,future,fut.var,maxent.path)
    } else {
      Biomod=fun.models.no.run(name,spdir)
    }
  }

  ##===================================================
  ## Plotting

  if (run.plots & !all(is.na(cell.pres))) {

    fun.plot(path,name,spdir,wcomp,p,zoom,enough,r.mar,e.map,Biomod[[1]],Biomod[[2]],fut.var,npix,environ,s,out.type)

  }

  ##===================================================
  ## Taxonomy

  if (run.taxo) {

    fun.taxo(path,name,spdir,spname,enough,npix)

  }

  return(enough)
}
