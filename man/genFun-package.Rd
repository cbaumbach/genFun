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
chr               Extract chromosome number from file names
chr2int           Convert chromosome string to integer
chunk             Extract chunk number from file names
combine_snptest   Summarize and combine snptest output by chromosome
extract_snps      Extract snps from imputed chunk files
find_genes        Map SNPs to genes
int2chr           Convert chromosome from integer to string
snptest           Run snptest gwas on multiple cores
unslash           Remove trailing slash from directory names
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

## combine_snptest
\dontrun{
#            OLD NAME                 NEW NAME
old2new <- c(rsid                   = "snp",
             chromosome             = "chr",
             position               = "pos",
                                      "imputed",
                                      "callrate",
             freq_alleleB           = "coding_freq",
             frequentist_add_beta_1 = "beta",
             frequentist_add_se_1   = "se",
             frequentist_add_pvalue = "pval")

# Custom `hook' function.
drop_snps_without_pvalue <- function(d)
{
    d[!is.na(d$pval), , drop = FALSE]
}

combine_snptest(indir    = "path/to/snptest_output_files",
                outdir   = "dir/where/output/should_go",
                ncore    = 6L,
                old2new  = old2new,
                select   = TRUE,
                template = "cohort_consortium_date_chr<CHROMOSOME>.txt.gz",
                gzip     = TRUE,
                hook     = drop_snps_without_pvalue)
}

## extract_snps
\dontrun{
extract_snps(
    snps      = c("rs123", "rs456"),
    indir     = "path/to/directory/with_chunk_files/",
    outdir    = "path/where/extracted/chunk/files/should/go",
    chunkmap  = c("path/to/chunkmap_part_1.txt", "path/to/chunkmap_part_2.txt"),
    idfile    = "path/to/file/with_all_person_ids_in_correct_order.txt",
    keep      = c("person_1", "person_2"),
    chunkmap_cols = c(1L, 3L, 4L),
    ncore     = 6L)
}

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
## Run snptest on men-only nonPAR chunks on chromosome X.
##
## Note: - The sample file must contain the same individuals as the
##         imputed chunk files.
##       - You have to specify "-assume_chromosome X".
##       - You have to craft your regular expression for `chr_chunk'
##         because the nonPAR chunk files don't follow the default
##         naming convention where chromosome and chunk are specified as
##         in "chr3_7".
##       - You have to specify a custom `pattern' because you need to
##         restrict the set of chunk files to those that contain only
##         male subjects.
snptest(
    indir          = "path/to/directory/with/gz-compressed/input/files",
    sample_file    = "path/to/MEN_ONLY_phenotype_file",
    outdir         = "path/to/directory/for/snptest/output/files",
    pheno          = "blood_sugar",
    covs           = c("age", "bmi"),
    exclusion_file = "path/to/file/with/exclusion/ids",
    add_args       = c("-missing_code NA",
                       "-frequentist 1",
                       "-method expected",
                       "-hwe",
                       "-printids",
                       "-use_raw_covariates",
                       "-use_raw_phenotypes",
                       "-assume_chromosome XY"),
    ncore          = 9L,
    pattern        = "males.*\\\\.gz$",
    chr_chunk      = ".*chr([^_]+)_(nonPAR_\\\\d+)",
    executable     = "/opt/snptest_v2.5_linux_x86_64_static/snptest_v2.5")
}

## unslash
unslash(c("path/to/dir_with_slash/", "path/to/dir_without_slash"))
}
