\name{extract_snps}
\alias{extract_snps}
\docType{package}
\title{
  Extract snps from imputed chunk files
}
\description{
  Extract snps from imputed chunk files
}
\usage{
extract_snps(snps, indir, chunkmap, chunkmap_cols = 1:3, keep = NULL,
             idfile = NULL, pattern = "\\\\.gz$", ncore = 1L,
             chr_chunk = ".*chr([^_]+)_(\\\\d+)")
}
\arguments{
  \item{snps}{
    Character vector with identifiers of snps to be extracted.
  }
  \item{indir}{
    Directory with gz-compressed chunk files.
  }
  \item{chunkmap}{
    Paths to files mapping snps to chromosomes and chunks.
  }
  \item{chunkmap_cols}{
    Length-3 integer vector naming the positions of the columns in the
    `chunkmap' files corresponding to snp identifier, chromosome, and
    chunk number.
  }
  \item{keep}{
    Character vector with identifiers of individuals for which to
    extract snps.
  }
  \item{idfile}{
    Path to file listing individual identifiers in the order in which
    their snp data appears in the imputed chunk files.
  }
  \item{pattern}{
    Regex pattern used to match input chunk files in `indir'.
  }
  \item{ncore}{
    Number of cores to use in parallel.
  }
  \item{chr_chunk}{
    Extended regular expression with two parenthesized subexpressions
    matching chromosome and chunk number in input chunk file names.
  }
}
\details{
  By default, the first 3 columns of the `chunkmap' files are assumed to
  correspond to the snp identifier, the chromosome, and the chunk number
  (in that order).  If your `chunkmap' files use different columns, you
  must specify the corresponding columns in `chunkmap_cols'.  Columns
  other than the ones named in `chunkmap_cols' are ignored.

  We assume the a `chunkmap' file has a header if all fields in the
  first row are of type "character".  Therefore, headers of your
  `chunkmap' files must not contain any numbers in the column names.

  The `idfile' must not have a header and must contain one identifier
  per line.
}
\value{
  Returns a data frame containing data for all snps in `snps' that could be
  extracted from the imputed chunk files in `indir'.  In case some snps
  could not be found, a warning is issued.
}
\examples{
library(genFun)

\dontrun{
extract_snps(
    snps      = c("rs123", "rs456"),
    indir     = "path/to/directory/with_chunk_files/",
    chunkmap  = c("path/to/chunkmap_part_1.txt", "path/to/chunkmap_part_2.txt"),
    chunkmap_cols = c(1L, 3L, 4L),
    keep      = c("person_1", "person_2"),
    idfile    = "path/to/file/with_all_person_ids_in_correct_order.txt",
    pattern   = "\\\\.gz$",
    ncore     = 1L,
    chr_chunk = ".*chr([^_]+)_(\\\\d+)"
}
}