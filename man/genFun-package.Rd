\name{genFun-package}
\alias{genFun-package}
\alias{genFun}
\docType{package}
\title{
Functions for dealing with genetic data
}
\description{
This is a collection of functions which, if not for this package, I
would have to rewrite or copy/paste for almost every R script dealing
with genetic data that I write.
}
\details{
\tabular{ll}{
Package: \tab genFun\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-04\cr
License: \tab GPL (>= 2)\cr
Packaged: \tab 2015-02-13 20:23:00 UTC; cb\cr
Built: \tab R 3.1.1; ; 2015-02-13 20:24:05 UTC; unix\cr
}

Index:
\preformatted{
find_genes       Map SNPs to genes
}
}
\author{
Clemens Baumbach

Maintainer: Clemens Baumbach <clemens.baumbach@helmholtz-muenchen.de>
}
\references{
none
}
\keyword{ manip }
\keyword{ attribute }
\keyword{ list }
\keyword{ utilities }
\keyword{ package }
\seealso{
  \code{\link[miscFun]{match_intervals}}
}
\examples{
library(genFun)

## find_genes
d <- read.table(textConnection("\
snp chr pos
rs1 1 1
rs2 1 5
rs3 2 3
rs4 3 4
", "r"), header = TRUE, stringsAsFactors = FALSE)

genes <- read.table(textConnection("\
id chr start end
g1 1 1 3
g2 1 1 6
g3 2 1 4
g4 4 1 5
", "r"), header = TRUE, stringsAsFactors = FALSE)

find_genes(d, genes)
}
