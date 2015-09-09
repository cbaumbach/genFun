\name{region_plot}
\alias{region_plot}
\docType{package}
\title{
  Plot snps with gene annotations
}
\description{
  Plot snps with gene annotations
}
\usage{
region_plot(file, x, y, start, end, label, width, ymax = NULL,
    transform = NULL, col = NULL, main = NULL)
}
\arguments{
  \item{file}{
    Name of pdf output file.
  }
  \item{x, y}{
    x and y coordinates of points to plot.
  }
  \item{start, end}{
    x coordinates of start and end positions of regions.
  }
  \item{label}{
    Labels for regions.
  }
  \item{width}{
    Width in inches of pdf output.
  }
  \item{ymax}{
    Value to replace infinite y-coordinates with.
  }
  \item{transform}{
    Function to be applied to y-coordinates before plotting.
  }
  \item{col}{
    Color for points (defaults to black).
  }
  \item{main}{
    String used as plot title.
  }
}
\value{
  None.
}
