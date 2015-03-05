\name{find_genes}
\alias{find_genes}
\title{
  Map SNPs to genes
}
\description{
  Given a data frame with SNPs (CpGs, etc.) and a data frame with genes
  return the first data frame with an additional column that contains
  for every SNP a comma-separated list with the identifiers of genes
  that contain the given SNP.
}
\usage{
find_genes(d, genes, chr1 = "chr", pos = "pos", out = "genes",
           chr2 = chr1, start = "start", end = "end", id = "id")
}
\arguments{
  \item{d}{
    Data frame with SNPs.
  }
  \item{genes}{
    Data frame with genes.
  }
  \item{chr1}{
    Name of column in `d' containing chromosome number.
  }
  \item{pos}{
    Name of column in `d' containing basepair position.
  }
  \item{out}{
    Name of column in output data frame that contains the
    comma-separated lists of genes.
  }
  \item{chr2}{
    Name of column in `genes' containing chromosome number.
  }
  \item{start}{
    Name of column in `genes' containing gene start positions.
  }
  \item{end}{
    Name of column in `genes' containing gene end positions.
  }
  \item{id}{
    Name of column in `genes' containing gene identifiers that will
    appear in the comma-separated lists of the new column `out' in the
    output data frame.
  }
}
\details{
  Both `d' and `genes' can contain other columns beside `chr1', `pos',
  `chr2', `start', `end', and `id'.  The order of SNPs and genes within
  the data frames is of no importance.  If a SNP does not belong to any
  of the genes, it will get a value of NA in the `out' column of the
  resulting data frame.  NAs in relevant columns of the SNP data frame
  lead to NAs in `out' column in output.  NAs in genes data frame are
  dropped before SNPs are matched to genes.
}
\value{
  Returns a data frame like `d' but with an additional column called
  `out' that contains for every SNP in `d' a comma-separated list of
  identifiers of genes from the `id' column in `genes' that contain the
  given SNP.
}
\seealso{
  \code{\link[miscFun]{match_intervals}}
}
\examples{
library(genFun)

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
