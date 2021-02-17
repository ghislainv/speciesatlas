# library(devtools)
# devtools::install_github("Warkorentin/speciesatlas")
# devtools::install("~/ANALYSE/Atlas_Corentin/speciesatlas",upgrade="never")
devtools::install("/home/corentinknoploch/Documents/Nitidae/Speciatlas/SpeciesAtlasPackage/",upgrade="never")

setwd("/home/corentinknoploch/Documents/Nitidae/Speciatlas/Lemurs/")
#
listPackages=data.frame()
listPackages <- c("gdalUtils","rgdal","raster","biomod2","ggplot2","knitr","xtable","magick","readbitmap","curl","htm2txt","dplyr","taxize","foreach","doParallel","parallel","bookdown","sp","grid")
lapply(listPackages, library, character.only = TRUE)

library(speciesatlas)
library(raster)

run.models <- T
run.plots <- T
run.taxo <- F
run.map <- F

#taxon.names <- c("lemurs","baobabs","birds")  "reptiles","geckos","trees"
taxon.names <- c("lemurs")
df.orig <- vector("list", length = length(taxon.names))
#taxon <- read.csv(file=paste0("data/",taxon.names[sp],"_prototype.txt"),header=TRUE,sep=";")

for (sp in 1:length(taxon.names)){
  taxon <- read.csv(file=paste0("data/",taxon.names[sp],".txt"),header=TRUE,sep=";")
  # taxon <- taxon[taxon$species %in% unlist(lapply(list.files("~/ANALYSE/Atlas_Corentin/atlas-run/figures/birds"), FUN=function(x) {gsub("\\.", " ",x)})),] #Randomly selecting 50 species
  # df.orig[[sp]]$Species <- droplevels(as.factor(taxon$species))
  df.orig[[sp]]$Species <- taxon$ScientificName
  df.orig[[sp]]$Lat <- taxon$decimallatitude
  df.orig[[sp]]$Long <- taxon$decimallongitude
  df.orig[[sp]]$Taxo.group <- taxon.names[sp]
}

##=======================
## Environmental data
##=======================

model.var <-  c("temp","prec","tseas","pseas","cwd","foret")

environ <- stack(stack("data/rastersNew/current.tif"),stack("data/rastersNew/environ.tif"))
names(environ) <- c(paste0("current",1:36),"temp",paste0("current",38:39),"tseas",paste0("current",41:45),
                    "meanWQ","meanCQ","prec","pWM","bioclim2","pseas",
                    paste0("current",52:67),"pet","cwd","ndm","alt",paste0("environ",2:7),"watershed","foret")
environ <- stack(environ) # Transform back from RasterBrick to rastertack

##=======================
## Future data
##=======================

fut.var <- list(c("cc","he","gs"),c("00"),c("2080"))
future <- vector("list")
i.mod <- 1

for (j in 1:length(fut.var[[2]])) { #rcp
  for (l in 1:length(fut.var[[3]])) { #yr
    future[[i.mod]] <- vector("list")
    environFut <- stack(paste0("data/rastersNew/environ_",fut.var[[3]][l],".tif"))
    if(fut.var[[2]][j]=="00"){
      for (mc in 1:length(fut.var[[1]])) {
        future[[i.mod]][[mc]]<- stack(dropLayer(environ,(nlayers(environ)-8):nlayers(environ)),environFut)
        names(future[[i.mod]][[mc]])<-c(paste0("current",1:36),"temp",paste0("current",38:39),"tseas",paste0("current",41:45),
            "meanWQ","meanCQ","prec","pWM","bioclim2","pseas",
            paste0("current",52:67),"pet","cwd","ndm","alt",paste0("environ",2:7),"watershed","foret")
        future[[i.mod]][[mc]]<-stack(future[[i.mod]][[mc]])
      }
    }else{
      for (mc in 1:length(fut.var[[1]])) { #mod
        ## Load climatic data
        future[[i.mod]][[mc]] <- stack(paste0("data/rastersNew/",fut.var[[1]][mc],"_",fut.var[[2]][j],"_",fut.var[[3]][l],".tif"),environFut)
        names(future[[i.mod]][[mc]]) <-c(paste0("current",1:36),"temp",paste0("current",38:39),"tseas",paste0("current",41:45),
                                         "meanWQ","meanCQ","prec","pWM","bioclim2","pseas",
                                         paste0("current",52:67),"pet","cwd","ndm","alt",paste0("environ",2:7),"watershed","foret")
        future[[i.mod]][[mc]] <- stack(future[[i.mod]][[mc]])
      } 
    }
    i.mod <- i.mod+1
  }
}

n.core <- 7
out.type="pdf"
title.book <- " Prototype of the Biodiversity Altas of Madagascar and its vulnerability to climate change"
author.book <- "G. Vieilledent $^a$, C.Knoploch $^b$, C. Grinand $^c$, F. Montfort $^c$, M. Nourtier $^c$"

fun.main(df.orig,run.models,run.plots,run.taxo,run.map,model.var,environ,future,fut.var,n.core,out.type,title.book,author.book)
