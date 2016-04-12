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
    in_terms_of_allele = c("B", "A"), transpose = FALSE)
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
  \item{transpose}{
    Should the resulting matrix should have individuals as rows and
    dosages as columns?
  }
}
\details{
Rows of the input matrix correspond to snps, consecutive triples of
columns to individuals.  Genotype probabilities are assumed to appear
in the order AA, AB, BB.  If `transpose' is FALSE, the resulting
matrix will have one row per snp and one column per individual.  If
`transpose' is TRUE, the resulting matrix will have individuals in the
rows and snps in the columns.  Both snps and individuals will appear
in the same order as in the input matrix.
}
\value{
Matrix with dosages.
}
\examples{
convert_to_dosages(rbind(.1,.2,.7))
}
