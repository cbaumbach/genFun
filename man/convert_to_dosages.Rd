\name{convert_to_dosages}
\alias{convert_to_dosages}
\docType{package}
\title{
Convert matrix of genotype probabilities to dosages
}
\description{
Convert matrix of genotype probabilities to dosages
}
\usage{
convert_to_dosages(aMatrix,
    model = c("additive", "recessive", "dominant"),
    in_terms_of_allele = c("B", "A"))
}
\arguments{
  \item{aMatrix}{
    Numeric matrix of genotype probabilities.
  }
  \item{model}{
    Type of genetic model to be used.
  }
  \item{in_terms_of_allele}{
    Allele in terms of which dosages should be computed.
  }
}
\details{
Genotype probabilities are assumed to appear in the order AA, AB, BB.
}
\value{
Matrix with dosages.
}
\examples{
convert_to_dosages(rbind(.1,.2,.7))
}
