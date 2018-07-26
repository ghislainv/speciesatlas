#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurelien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

fun.main <- function(df.orig,run.models=TRUE,run.plots=TRUE,run.taxo=TRUE,run.map=TRUE,model.var,environ,future,fut.var,maxent.path,n.core=(detectCores()-1),out.type="html",title.book="Title",author.book="Author"){

  # =======================
  # API Keys
  # =======================

  ## API key for taxize (check yours are in an Renviron.site file in the wd)
  ## See: https://ropensci.github.io/taxize-book/authentication.html
  eol.key <- getkey(service="eol")
  tropicos.key <- getkey(service="tropicos")
  iucn.key <- getkey(service="iucn")
  ncbi.key <- getkey(service="entrez")
  check_APIkeys <- c(eol.key,tropicos.key,iucn.key,ncbi.key)

  SP <- fun.data(df.orig)
  df.sp <- SP[[1]]
  sp.names <- SP[[2]]
  sp.dir <- SP[[3]]
  taxon.sp <- SP[[4]]
  taxon.names <- SP[[5]]

  ##=======================
  ## Environmental data
  s <- stack(environ[[model.var]])

  ##=======================
  ## Parallel computations
  ##=======================

  ## For MAXENT.Phillips with JAVA to work on RStudio server
  Sys.unsetenv("DISPLAY")
  ## Package names for parallel computations
  pkg.names.clust <- c("rgdal","raster","biomod2","ggplot2","knitr","grid",
                       "xtable","magick","readbitmap","curl","htm2txt","dplyr","taxize",
                       "foreach","doParallel","parallel","bookdown","sp")

  ## Make a cluster with all possible cores
  if(n.core>detectCores()){n.core <- detectCores()-1}
  clust <- makeCluster(n.core)
  ## Register the number of parallel workers (here all CPUs)
  registerDoParallel(clust)
  ## Return number of parallel workers
  getDoParWorkers()
  enough <- vector() ## Retain if there is enough observations
  enough <- foreach(i=1:length(taxon.sp),.packages=pkg.names.clust) %dopar% fun.species(i,run.models,run.plots,run.taxo,model.var,environ,future,fut.var,maxent.path,out.type,taxon.sp,taxon.names,sp.dir,sp.names,df.sp,s)
  ## Stop the cluster
  stopCluster(clust)
  setwd("..")

  if(run.map){
    fun.map(sp.dir,extent(environ),enough,taxon.names,taxon.sp,fut.var)
  }

  if((out.type=="html")||(out.type=="both")){
    fun.book(sp.dir,taxon.names,taxon.sp,enough)
  }
  if((out.type=="pdf")||(out.type=="both")){
    fun.pdf(sp.dir,taxon.names,taxon.sp,enough)
  }

}
