\name{chr2int}
\alias{chr2int}
\docType{package}
\title{
  Convert chromosome string to integer
}
\description{
  Convert chromosome string to integer.
}
\usage{
chr2int(x, prefix = NULL)
}
\arguments{
  \item{x}{
    Character vector of chromosome strings to be converted to integer.
  }
  \item{prefix}{
    Prefix to be removed from chromosome strings before conversion.
  }
}
\value{
  Returns chromosome numbers associated with chromosome strings in x as
  integers.
}
\examples{
library(genFun)

chr2int(c("1", "MT", "X", "Y"))
chr2int(c("chr1", "chrMT", "chrX", "chrY"), prefix = "chr")
}
