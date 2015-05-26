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
snptest(indir, sample_file, outdir, pheno, covs = NULL,
        exclusion_file = NULL, add_args = NULL, ncore = 1L,
        pattern = "\\\\.gz$", chr_chunk = ".*chr([^_]+)_(\\\\d+)",
        executable = "snptest")
}
\arguments{
  \item{indir}{
    Directory with gz-compressed chunk files.
  }
  \item{sample_file}{
    Phenotype file.
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
  \item{exclusion_file}{
    File with ids of individuals to be excluded from analysis.
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
  \item{chr_chunk}{
    Extended regular expression with two parenthesized subexpressions
    matching chromosome and chunk number in input chunk file names.
  }
  \item{executable}{
    Path to snptest executable.
  }
}
\details{
  Function parameters without default values are mandatory.  All other
  parameters are optional.

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

  If your chunk files follow a naming scheme where chromosome and chunk
  number cannot be matched by the default regular expression for
  `chr_chunk', then you will have to specify your own regular expression
  (see \link[base]{regex}).  This (extended) regular expression must
  contain two parenthesized subexpressions.  The first must match the
  chromosome, the second the chunk part of your chunk file names.
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
}
