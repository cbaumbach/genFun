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
    transform = NULL, col = NULL, main = NULL, upper_hook = NULL,
    lower_hook = NULL, point_cex = .3)
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
  \item{upper_hook}{
    Function that will be called right after the upper part of the plot
    is created.  Its arguments correspond to `x', `y', `start', `end',
    and `label'.
  }
  \item{lower_hook}{
    Same as `upper_hook' except that the function is called right after
    the lower part of the plot is created.
  }
  \item{point_cex}{
    Factor by which points in plot are magnified.
  }
}
\value{
  None.
}
