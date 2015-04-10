\name{unslash}
\alias{unslash}
\docType{package}
\title{
  Remove trailing slash from directory names
}
\description{
  Remove trailing slash from directory names
}
\usage{
unslash(dirs)
}
\arguments{
  \item{dirs}{
    Character vector of directory paths.
  }
}
\value{
  Returns `dirs' with any trailing slashes removed.
}
\examples{
library(genFun)

unslash(c("path/to/dir_with_slash/", "path/to/dir_without_slash"))
}
