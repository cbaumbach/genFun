\name{snptest}
\alias{snptest}
\docType{package}
\title{
  Run snptest gwas on multiple cores
}
\description{
  Run snptest gwas on multiple cores
}
\usage{
snptest(indir, sample_file, exclusion_file, outdir, pheno, covs,
        add_args = NULL, ncore = 1L, pattern = "\\\\.gz$",
        chunk_pattern = NULL, executable = "snptest")
}
\arguments{
  \item{indir}{
    Directory with gz-compressed chunk files.
  }
  \item{sample_file}{
    Phenotype file.
  }
  \item{exclusion_file}{
    File with ids of individuals to be excluded from analysis.
  }
  \item{outdir}{
    Directory where snptest output files should go.
  }
  \item{pheno}{
    Name of phenotype variable as used in `sample_file'.
  }
  \item{covs}{
    Character vector with names of covariates as used in `sample_file'.
  }
  \item{add_args}{
    Additional command-line arguments to snptest.
  }
  \item{ncore}{
    Number of cores to use in parallel.
  }
  \item{pattern}{
    Regex pattern used to match input chunk files in `indir'.
  }
  \item{chunk_pattern}{
    Allows the user to specify a custom (extended) regular expression in
    which the chunk part of the input file names is matched by a
    parenthesized subexpression (e.g. "chr[^_]+_([[:alpha:]]+_\\\\d+)").
  }
  \item{executable}{
    Path to snptest executable.
  }
}
\details{
  Sometimes you do not want all chunk files in `indir' to be analyzed.
  A typical reason would be that different chunk files need to be run
  with different snptest options.  If you are able to specify a regular
  expression that matches the set of chunk files that you do want to
  analyze, you can supply this regular expression as the `pattern'
  argument to \code{snptest}.  This will cause snptest to be run only
  with files that are located in `indir' _and_ whose path---which
  includes `indir' as a prefix (!)---matches the regular expression
  `pattern'.  Make sure that your pattern a) matches only chunk files
  and b) that all matched chunk files can be snptest'ed using the same
  set of snptest options.

  If your chunk files follow a naming scheme where the chunk number
  cannot be matched by the parenthesized subexpression of the following
  regular expression "chr[^_]+_(\\\\d+)", then you will have to specify
  your own regular expression (see \link[base]{regex}).  This (extended)
  regular expression must contain a parenthesized subexpression which
  matches the chunk part of your chunk file names.
}
\value{
  Returns a data frame with columns: `chr' (chromosome number), `chunk'
  (chunk number), `input' (gz-compressed input file), `output' (snptest
  output file), `log' (snptest log file), `done' (TRUE for successfully
  snptest'ed chunks, otherwise FALSE).

  After an entirely successful run of \code{snptest} all values
  in the `done' column will be "TRUE".  Otherwise columns with `done ==
  "FALSE"' should give you a hint where to look for possible problems.
  Also take a look at the warning, if any, that \code{snptest}
  emits while running.
}
\examples{
library(genFun)

\dontrun{
## Run snptest on men-only nonPAR chunks on chromosome X.
##
## Note: - The sample file must contain the same individuals as the
##         imputed chunk files.
##       - You have to specify "-assume_chromosome X".
##       - You have to use a custom `chunk_pattern' because the nonPAR
##         chunk files don't follow the default naming convention where
##         chromosome and chunk are specified as in "chr3_7".
##       - You have to specify a custom `pattern' because you need to
##         restrict the set of chunk files to those that contain only
##         male subjects.
snptest(
    indir          = "path/to/directory/with/gz-compressed/input/files",
    sample_file    = "path/to/MEN_ONLY_phenotype_file",
    exclusion_file = "path/to/file/with/exclusion/ids",
    outdir         = "path/to/directory/for/snptest/output/files",
    pheno          = "blood_sugar",
    covs           = c("age", "bmi"),
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
    chunk_pattern  = "chr[^_]+_(nonPAR_\\\\d+)",
    executable     = "/opt/snptest_v2.5_linux_x86_64_static/snptest_v2.5")
}
}
