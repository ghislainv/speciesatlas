#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurelien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

##=======================
## Render the book
##=======================

fun.book <- function(sp.dir,taxon.names,taxon.sp,enough){
  rmd.path <- system.file("Rmd", package="speciesatlas")

  file.copy(paste(rmd.path,"index.Rmd",sep = "/"),getwd())
  file.copy(paste(rmd.path,"atlas.sp.Rmkd",sep = "/"),getwd())
  file.copy(paste(rmd.path,"atlas.sp.not.enough.Rmkd",sep = "/"),getwd())
  file.copy(paste(rmd.path,"taxon.Rmkd",sep = "/"),getwd())
  file.copy(paste(rmd.path,"style.css",sep = "/"),getwd())
  file.copy(paste(rmd.path,"_bookdown.yml",sep = "/"),getwd())
  file.copy(paste(rmd.path,"doc_prefix.tex",sep = "/"),getwd())
  file.copy(paste(rmd.path,"header.tex",sep = "/"),getwd())
  file.copy(paste(rmd.path,"biblio.bib",sep = "/"),getwd())
  file.copy(paste(rmd.path,"jae.bst",sep = "/"),getwd())
  file.copy(paste(rmd.path,"journal-of-applied-ecology.csl",sep = "/"),getwd())
  file.copy(paste(rmd.path,"../img/head.png",sep = "/"),paste0(getwd(),"/figures"))
  file.copy(paste(rmd.path,"../img/alt.png",sep = "/"),paste0(getwd(),"/figures"))
  file.copy(paste(rmd.path,"../img/ca.png",sep = "/"),paste0(getwd(),"/figures"))
  file.copy(paste(rmd.path,"../img/niche.png",sep = "/"),paste0(getwd(),"/figures"))
  file.copy(paste(rmd.path,"../img/cafd.png",sep = "/"),paste0(getwd(),"/figures"))
  file.copy(paste(rmd.path,"../img/cazd.png",sep = "/"),paste0(getwd(),"/figures"))
  file.copy(paste(rmd.path,"../img/plotting.rda",sep = "/"),paste0(getwd(),"/figures"))

  wSp <- 1:length(taxon.sp)
  # html
  # Don't indicate output_format to take into account YAML options
  options(knitr.table.format="html")
  # Dynamic YAML options
  params <- list(title=title.book,author=author.book,date=format(Sys.time(), "%d %B, %Y"))
  bookdown::render_book("index.Rmd")


  file.remove("index.Rmd")
  file.remove("atlas.sp.Rmkd")
  file.remove("atlas.sp.not.enough.Rmkd")
  file.remove("taxon.Rmkd")
  file.remove("style.css")
  file.remove("_bookdown.yml")
  file.remove("doc_prefix.tex")
  file.remove("header.tex")
  file.remove("biblio.bib")
  file.remove("jae.bst")
  file.remove("journal-of-applied-ecology.csl")
  file.remove("figures/head.png")
  file.remove("figures/alt.png")
  file.remove("figures/ca.png")
  file.remove("figures/niche.png")
  file.remove("figures/cafd.png")
  file.remove("figures/cazd.png")
  file.remove("figures/plotting.rda")
}
