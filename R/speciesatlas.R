#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurélien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

fun.main <- function(df.orig,run.models,run.plots,run.taxo,environ,future,fut.var,maxent.path,n.core,title.book,author.book){

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
  ## Parallel computations
  ##=======================

  ## For MAXENT.Phillips with JAVA to work on RStudio server
  Sys.unsetenv("DISPLAY")
  ## Use bytecode compilation
  fun.species.bytes <- cmpfun(fun.species)
  ## Package names for parallel computations
  pkg.names.clust <- c("rgdal","raster","biomod2","ggplot2","Reol","fields","knitr",
                       "xtable","magick","readbitmap","curl","htm2txt","dplyr","taxize",
                       "foreach","doParallel","parallel","compiler","bookdown","kableExtra","sp")

  ## Make a cluster with all possible cores
  if(n.core>detectCores()){n.core <- detectCores()-1}
  clust <- makeCluster(n.core)
  ## Register the number of parallel workers (here all CPUs)
  registerDoParallel(clust)
  ## Return number of parallel workers
  getDoParWorkers()
  enough <- vector() ## Retain if there is enough observations
  enough <- foreach(i=1:length(taxon.sp),.packages=pkg.names.clust) %dopar% fun.species.bytes(i,run.models,run.plots,run.taxo,environ,future,fut.var,maxent.path,taxon.sp,taxon.names,sp.dir,sp.names,df.sp)
  ## Stop the cluster
  stopCluster(clust)
  setwd("..")

  fun.book(sp.dir,taxon.names,taxon.sp,enough)
}
