\name{read.obo}
\alias{read.obo}
\docType{package}
\title{
  Read gene ontology obo 1.2 file
}
\description{
  Read gene ontology obo 1.2 file
}
\usage{
read.obo(filename, colClasses = NULL)
}
\arguments{
  \item{filename}{
    Path to obo file.
  }
  \item{colClasses}{
    Character vector as used by \code{read.table}.
  }
}
\value{
  Returns a data frame with one go term per row and the term-specific
  values of selected obo tags as columns.
}
\seealso{
  \code{\link[utils]{read.table}}
}
\examples{
library(genFun)

\dontrun{
cat("\
[Term]
id: GO:0000000
name: foobar
alt_id: GO:7654321
is_a: GO:1111111
is_a: GO:2222222
def: \"An imaginary go term.\"

[Term]
id: GO:3333333
name: bazqux
is_a: GO:4444444
def: \"Another imaginary go term.\"
", file = "foo.obo")
d <- read.obo("foo.obo")
}
}
