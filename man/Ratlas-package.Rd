\name{speciesatlas-package}
\alias{speciesatlas-package}
\alias{speciesatlas}
\docType{package}

\title{Package speciesatlas}

\description{\packageDescription{speciesatlas}}

\details{
To use this package, you must have :

- a list of dataFrames containing your data. Each dataFrame stands for a taxon and must contains at least these columns with the exact same names : "Species" (name of the specie), "Long" (Longitude), "Lat" (Latitude) and "Taxo.group" (name of the taxon).

- a RasterStack containing environmental data about the area you want to study. It has to contain an "alt" column, and the other columns are different parameters. You decide which ones you want to consider runing the models.

- a list of RasterStacks with environmental datas about the future, according to several climate scenarios. For example if you have 2 years, 2 rcp scenarios and 3 climate models, your object will be a list of 4 list with 3 rasters in each.

- your API Keys in an .Renviron file at the root of your computer. More informations at \url{https://ropensci.github.io/taxize-book/authentication.html}.

- a directory called "maxent" containing the file "maxent.jar"
}

\author{
\packageAuthor{speciesatlas}

Maintainer: \packageMaintainer{speciesatlas}}

\references{
Muniz-Tagliari M., J.-M. Leong Pock-Tsy, C. Cornu and P. Danthu and G. Vieilledent. Vulnerability of the seven baobab species in Madagascar to climate change. in prep.

Vieilledent G., M. Muniz-Tagliari, C. Grinand, F. Montfort. Atlas of the biodiversity of Madagascar: present species distribution and species vulnerability to climate change. in prep.
}

\keyword{speciesatlas}
\seealso{
\url{https://bioscenemada.cirad.fr/}

\url{https://github.com/ghislainv/atlas}
}
\examples{
fun.main(df.orig,run.models,run.plots,run.taxo,model.var,environ,future,fut.var,maxent.path,n.core,title.book,author.book)
}
