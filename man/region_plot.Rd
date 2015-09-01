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
    transform = NULL, col = NULL, main = NULL, special = NULL,
    special_col = NULL)
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
    Color used to mark non-special points.
  }
  \item{main}{
    String used as plot title.
  }
  \item{special}{
    Logical vector of same length as `x' and `y' indicating which points
    should be drawn using `special_col'.
  }
  \item{special_col}{
    Color used for points that are `special'.
  }
}
\value{
  None.
}
