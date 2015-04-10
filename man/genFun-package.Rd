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
chr              Extract chromosome number from file names
chr2int          Convert chromosome string to integer
chunk            Extract chunk number from file names
find_genes       Map SNPs to genes
int2chr          Convert chromosome from integer to string
snptest          Run snptest gwas on multiple cores
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

## chr
chr(c("chr5_path/to/the_chr1_3.txt", "chr5_path/to/the_chr2_4.txt"))

## chr2int
chr2int(c("1", "MT", "X", "Y"))
chr2int(c("chr1", "chrMT", "chrX", "chrY"), prefix = "chr")

## chunk
chunk(c("chr5_path/to/the_chr1_3.txt", "chr5_path/to/the_chr2_4.txt"))

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

## int2chr
int2chr(c(1:3, 23:25))

## snptest
\dontrun{
snptest(
    indir          = "path/to/directory/with/gz-compressed/input/files",
    sample_file    = "path/to/phenotype/file",
    exclusion_file = "path/to/file/with/exclusion/ids",
    outdir         = "path/to/directory/for/snptest/output/files",
    pheno          = "blood_sugar",
    covs           = c("age", "sex", "bmi"),
    add_args       = c("-missing_code NA",
                       "-frequentist 1",
                       "-method expected",
                       "-hwe",
                       "-printids",
                       "-use_raw_covariates",
                       "-use_raw_phenotypes"),
    ncore          = 9L,
    executable     = "/opt/snptest_v2.5_linux_x86_64_static/snptest_v2.5")
}

}
