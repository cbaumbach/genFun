\name{count_per_gene}
\alias{count_per_gene}
\docType{package}
\title{
Count number of points per gene
}
\description{
Count number of points per gene
}
\usage{
count_per_gene(pos, chr, gene_start, gene_end, gene_chr)
}
\arguments{
  \item{pos}{
    positions
  }
  \item{chr}{
    chromosome numbers of positions
  }
  \item{gene_start}{
    left endpoints of gene intervals
  }
  \item{gene_end}{
    right endpoints of gene intervals
  }
  \item{gene_chr}{
    chromosome numbers of gene intervals
  }
}
\value{
An integer vector of the same length as `gene_start` containing the
number of points falling into each of the gene intervals.
}
\examples{
stopifnot(count_per_gene(0, 0, 0, 0, 0) == 1)
}
