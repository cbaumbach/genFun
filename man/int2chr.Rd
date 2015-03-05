\name{int2chr}
\alias{int2chr}
\docType{package}
\title{
  Convert chromosome from integer to string
}
\description{
  Convert chromosome from integer to string.
}
\usage{
int2chr(x)
}
\arguments{
  \item{x}{
    Integer vector with chromosomes to be converted to strings.
  }
}
\value{
  Returns chromosomes associated with integer values in x as character
  vector.
}
\examples{
library(genFun)

int2chr(c(1:3, 23:25))
}
