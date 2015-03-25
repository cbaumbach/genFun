\name{dot_attributes}
\alias{dot_attributes}
\alias{print.dot_attributes}
\docType{package}
\title{
  Create and print dot format node or edge attributes
}
\description{
  Create and print dot format node or edge attributes
}
\usage{
dot_attributes(...)
\method{print}{dot_attributes}(x, ...)
}
\arguments{
  \item{x}{
    Object of type `dot_attributes'.
  }
  \item{...}{
    List of `key = value' pairs where the key corresponds to a dot
    format node or edge attribute.  Ignored for \code{print}.
  }
}
\value{
  \code{dot_attributes} returns a new object of type `dot_attributes'.
  \code{print} prints the list of dot attributes to standard output in
  the dot format.
}
\examples{
library(genFun)

print(dot_attributes(style = "dotted", color = "green"))
}
