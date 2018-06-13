#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurélien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

##=======================
## Render the book
##=======================

fun.book <- function(sp.dir,taxon.names,taxon.sp,enough){
  rmd_path <- system.file("Rmd", package="speciesatlas")

  file.copy(paste(rmd_path,"index.Rmd",sep = "/"),getwd())
  file.copy(paste(rmd_path,"atlas.sp.Rmkd",sep = "/"),getwd())
  file.copy(paste(rmd_path,"atlas.sp.not.enough.Rmkd",sep = "/"),getwd())
  file.copy(paste(rmd_path,"taxon.Rmkd",sep = "/"),getwd())
  file.copy(paste(rmd_path,"style.css",sep = "/"),getwd())
  file.copy(paste(rmd_path,"_bookdown.yml",sep = "/"),getwd())
  file.copy(paste(rmd_path,"doc_prefix.tex",sep = "/"),getwd())
  file.copy(paste(rmd_path,"header.tex",sep = "/"),getwd())
  file.copy(paste(rmd_path,"biblio.bib",sep = "/"),getwd())
  file.copy(paste(rmd_path,"jae.bst",sep = "/"),getwd())
  file.copy(paste(rmd_path,"journal-of-applied-ecology.csl",sep = "/"),getwd())

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
}
