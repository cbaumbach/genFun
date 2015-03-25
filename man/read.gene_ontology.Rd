\name{read.gene_ontology}
\alias{read.gene_ontology}
\docType{package}
\title{
  Create tree structure from gene ontology file
}
\description{
  Create tree structure from gene ontology file
}
\usage{
read.gene_ontology(filename, ontology = c(
                             "biological_process",
                             "cellular_component",
                             "molecular_function"))
}
\arguments{
  \item{filename}{
    Path to obo file with gene ontology.
  }
  \item{ontology}{
    Character vector of ontologies for which to create tree objects.
  }
}
\value{
  Returns a list with one tree object per `ontology'.
}
\seealso{
  \code{\link{read.obo}}
}
