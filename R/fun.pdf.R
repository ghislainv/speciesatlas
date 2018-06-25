#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurelien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

##=======================
## Render the pdf
##=======================

fun.pdf <- function(sp.dir,taxon.names,taxon.sp,enough){

  wSp <- 1:length(taxon.sp)
  ## Set knitr chunk default options
  opts_chunk$set(echo=FALSE, cache=FALSE,
                 results="hide", warning=FALSE,
                 message=FALSE, highlight=TRUE,
                 fig.show="hide", size="small",
                 tidy=FALSE)
  knitr::opts_knit$set(root.dir=normalizePath(getwd()))
  ## Knit
  user.path <- getwd()
  rnw.path <- system.file("Rnw", package="speciesatlas")
  dir.create("Pdf",recursive=TRUE,showWarnings=FALSE)
  knit2pdf(paste0(rnw.path,"/atlas.Rnw"),output="Pdf/atlas.tex")

}
