\name{chr}
\alias{chr}
\alias{chunk}
\docType{package}
\title{
  Extract chromosome or chunk from file names
}
\description{
  Extract chromosome or chunk from file names
}
\usage{
chr(xs, pattern = ".*chr([^_]+)")
chunk(xs, pattern = ".*chr[^_]+_(\\\\d+)")
}
\arguments{
  \item{xs}{
    File names.
  }
  \item{pattern}{
    Extended regular expression pattern (see \code{\link[base]{regex}})
    used to extract the substrings corresponding to chromosome or chunk
    number from `xs'.
  }
}
\value{
  \code{chr} returns an integer vector of chromosome numbers unless the
  conversion from character to integer would lead to missing values, in
  which case a character vector is returned.

  \code{chunk} returns an integer vector of chunk numbers.
}
\examples{
library(genFun)

chr(c("chr5_path/to/the_chr1_3.txt", "chr5_path/to/the_chr2_4.txt"))
chunk(c("chr5_path/to/the_chr1_3.txt", "chr5_path/to/the_chr2_4.txt"))
}
