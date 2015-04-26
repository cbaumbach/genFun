\name{combine_snptest}
\alias{combine_snptest}
\docType{package}
\title{
  Summarize and combine \code{snptest} output by chromosome
}
\description{
  \code{combine_snptest} proceeds by chromosome: it summarizes and
  combines all snptest output files for a chromosome and writes the
  resulting data frame to a file in `outdir'.  The chromosome-specific
  output file is named according to `template' with the "<CHROMOSOME>"
  placeholder replaced by the actual chromosome.
}
\usage{
combine_snptest(indir, outdir, ncore = 1L, old2new = NULL,
                select = !is.null(old2new), pattern = "^chr.*\\\\.txt$",
                template = "chr<CHROMOSOME>.txt",
                chr_chunk = ".*chr([^_]+)_(\\\\d+)", gzip = FALSE,
                hook = NULL)
}
\arguments{
  \item{indir}{
    Directory with snptest output files.
  }
  \item{outdir}{
    Directory where summarized and combined snptest output should go.
  }
  \item{ncore}{
    Number of cores to use in parallel.
  }
  \item{old2new}{
    Character vector used for renaming, selecting, and reordering
    columns from snptest output (see Details).
  }
  \item{select}{
    If TRUE, `old2new' will be used to select and reorder a subset of
    columns for snptest output.
  }
  \item{pattern}{
    Regex pattern used to match snptest output files in `indir'.
  }
  \item{template}{
    String serving as template for naming output files.  Must contain
    the substring "<CHROMOSOME>" which will be replaced by the actual
    chromosome.
  }
  \item{chr_chunk}{
    Extended regular expression with two parenthesized subexpressions
    matching chromosome and chunk number in snptest output files in
    `indir'.
  }
  \item{gzip}{
    If TRUE chromosome-specific output files will be compressed before
    being written to disk.
  }
  \item{hook}{
    Function accepting and returning a data frame (see Details).
  }
}
\details{
  The `old2new' argument, if non-NULL, serves up to three purposes.
  \enumerate{
    \item Named elements of `old2new' will be used to rename columns
          from the snptest output files.  An element of `old2new' with
          name `old_name' and value `new_name' will rename the
          `old_name' column in a snptest output file to `new_name'.
    \item The values of `old2new' determine which columns of the snptest
          output files will be included in the returned data frame.
    \item The order of the value in `old2new' determines the order of
          columns in the returned data frame.
  }

  Note that elements of `old2new' that are only used to select and
  reorder columns need not be named.  In other words, `old2new' may
  contain a mixture of named and unname elements.  In fact, if the
  `names' attribute of `old2new' is NULL and `select = TRUE', `old2new'
  will not be used for renaming but only for selecting and reordering
  columns.

  In addition to the columns that the snptest program creates in its
  output files there will be the following columns:
  \itemize{
    \item `freq_alleleB' gives the frequency of the allele in the
          `alleleB' column.  Since in snptest this is the effect allele,
          `freq_alleleB' is the effect allele frequency.  It is computed
          as (`all_AB + (2 * `all_BB')) / (2 * (`all_AA' + `all_AB' +
          `all_BB')).
    \item `callrate' gives the fraction of non-missing genotypes.  It is
          computed as 1 - (`all_NULL' / (`all_AA' + `all_AB' + `all_BB'
          + `all_NULL')).
    \item `imputed' is 1 for snps that were imputed and otherwise 0.
  }

  If `old2new' is NULL all columns from the snptest output files will be
  included in the chromosome-specific output files of
  \code{combine_snptest}, including the above-mentioned extra columns.

  Sometimes being able to rename, select, and reorder columns via
  `old2new' and `select' doesn't quite cut it.  If you need unlimited
  flexibility, you can supply a function via the `hook' argument.  The
  function must take and return a data frame.  It will be applied to a
  data frame representing the summarized contents of a single (!)
  snptest output file after columns have been renamed, selected, and
  reordered according to `old2new' and `select'.  Make sure that your
  `hook' function does indeed return a data frame.  In particular, use
  `drop = FALSE' when subsetting with \code{[}.
}
\value{
  None.
}
\examples{
library(genFun)

\dontrun{
# Define how to rename, select, and reorder columns.
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
}
