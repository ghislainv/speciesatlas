\name{fun.main}
\alias{fun.main}

\title{Main function of the speciesatlas package}

\description{
This function can be called to make a biodiversity atlas of a given area.
}

\usage{
fun.main(df.orig,run.models=TRUE,run.plots=TRUE,run.taxo=TRUE,model.var,environ,future,fut.var,maxent.path,n.core=(detectCores()-1),out.type="both",title.book="Title",author.book="Author",eol.key=getkey(service="eol"),tropicos.key=getkey(service="tropicos"),iucn.key=getkey(service="iucn"),ncbi.key=getkey(service="entrez"))
}

\arguments{
  \item{df.orig}{
List of dataFrames containing your data. Each dataFrame stands for a taxon and must contains at least these columns with the exact same names : "Species" (name of the specie), "Long" (Longitude), "Lat" (Latitude) and "Taxo.group" (name of the taxon).
}
  \item{run.models}{
Indicate if you want the models part to be run or not.
}
  \item{run.plots}{
Indicate if you want the plots part to be run or not.
}
  \item{run.taxo}{
Indicate if you want the taxonomy part to be run or not.
}
  \item{model.var}{
Contains the parameters you want to run the study with.
}
  \item{environ}{
RasterStack containing environmental data about the area you want to study. It has to contain an "alt" column, and the other columns are different parameters.
}
  \item{future}{
List of RasterStacks with environmental datas about the future, according to several climate scenarios.For example if you have 2 years, 2 rcp scenarios and 3 climate models, your object will be a list of 4 list with 3 rasters in each.

Exemple of future list :

future <- list("cc_45_2050","cc_45_2080","cc_85_2050","cc_85_2080","he_45_2050","he_45_2080","he_85_2050","he_85_2080")
}
  \item{fut.var}{
List of the future variables you chose. Put the different climate models first, afterthat the different rcp scenarios and finally the different years.

Exemple of fut.var list :
fut.var <- list(c("cc","he"),c("45","85"),c("2050","2080"))
}
  \item{maxent.path}{
Directory where maxent.jar can be found.
}
  \item{n.core}{
Number of cores you want to be used by running the models
}
  \item{title.book}{
The title you want for your book.
}
  \item{author.book}{
Authors of your book.
}
  \item{eol.key}{
Your Encyclopedia of Life API Key. The keys can be stored in a .REnviron file, as described in the details section.
}
  \item{tropicos.key}{
Your Tropicos API Key. The keys can be stored in a .REnviron file, as described in the details section.
}
  \item{iucn.key}{
Your IUCN API Key. The keys can be stored in a .REnviron file, as described in the details section.
}
  \item{ncbi.key}{
Your NCBI API Key. The keys can be stored in a .REnviron file, as described in the details section.
}
}

\details{
You must have your API Keys in an .Renviron file at the root of your computer. More informations at \url{https://ropensci.github.io/taxize-book/authentication.html}.

Also a directory called "maxent" must contain the file "maxent.jar", the path to it is given in the arguments.
}

\value{
This function produces an output in the Bookdown format. It is created in the "Book" directory at the root.
}

\references{
Muniz-Tagliari M., J.-M. Leong Pock-Tsy, C. Cornu and P. Danthu and G. Vieilledent. Vulnerability of the seven baobab species in Madagascar to climate change. in prep.

Vieilledent G., M. Muniz-Tagliari, C. Grinand, F. Montfort. Atlas of the biodiversity of Madagascar: present species distribution and species vulnerability to climate change. in prep.
}

\author{
\packageAuthor{speciesatlas}
}

\seealso{
\url{https://bioscenemada.cirad.fr/}

\url{https://github.com/ghislainv/atlas}
}

\examples{

library(raster)

run.models <- FALSE
run.plots <- TRUE
run.taxo <- TRUE
df.orig <- list(read.csv(file="data/lemurs.txt",header=TRUE),read.csv(file="data/baobabs.txt",header=TRUE))

model.var <-  c("temp","tseas","prec","cwd")
environ <- stack("data/rasters/environ.tif")

fut.var <- list(c("cc","he","gs"),c("45","85"),c("2050","2080"))
future <- vector("list")
i.mod <- 1

for (j in 1:length(fut.var[[2]])) { #rcp
  for (l in 1:length(fut.var[[3]])) { #yr
    future[[i.mod]] <- vector("list")
    for (mc in 1:length(fut.var[[1]])) { #mod
      ## Load climatic data
      future[[i.mod]][[mc]] <- stack(paste0("data/rasters/",fut.var[[1]][mc],"_",fut.var[[2]][j],"_",fut.var[[3]][l],".tif"))
    }
    i.mod <- i.mod+1
  }
}

maxent.path <- paste0(getwd(),"/maxent")
n.core <- 3
title.book <- "Biodiversity Atlas of Madagascar"
author.book <- "G.Vieilledent, A.Colas"

eol.key <- "xxx"
tropicos.key <- "xxx"
iucn.key <- "xxx"
ncbi.key <- "xxx"

fun.main(df.orig,run.models,run.plots,run.taxo,environ,future,fut.var,maxent.path,n.core,title.book,author.book,eol.key,tropicos.key,iucn.key,ncbi.key)
}

\keyword{speciesatlas}
