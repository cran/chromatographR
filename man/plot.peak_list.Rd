% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_peaks.R
\name{plot.peak_list}
\alias{plot.peak_list}
\title{Plot fitted peak shapes.}
\usage{
\method{plot}{peak_list}(
  x,
  ...,
  chrom_list = NULL,
  index = 1,
  lambda = NULL,
  points = FALSE,
  ticks = FALSE,
  a = 0.5,
  color = NULL,
  cex.points = 0.5
)
}
\arguments{
\item{x}{Peak_list object. Output from the \code{get_peaks} function.}

\item{\dots}{Additional arguments to plot function.}

\item{chrom_list}{List of chromatograms (retention time x wavelength
matrices)}

\item{index}{Index or name of chromatogram to be plotted.}

\item{lambda}{Wavelength for plotting.}

\item{points}{Logical. If TRUE, plot peak maxima. Defaults to FALSE.}

\item{ticks}{Logical. If TRUE, mark beginning and end of each peak. Defaults
to FALSE.}

\item{a}{Alpha parameter controlling the transparency of fitted shapes.}

\item{color}{The color of the fitted shapes.}

\item{cex.points}{Size of points. Defaults to 0.5}
}
\value{
No return value, called for side effects.
}
\description{
Visually assess integration accuracy by plotting fitted peaks over trace.
}
\section{Side effects}{

Plots a chromatographic trace from the specified chromatogram (\code{chr})
at the specified wavelength (\code{lambda}) with fitted peak shapes from the
provided \code{peak_list} drawn underneath the curve.
}

\seealso{
\code{\link{get_peaks}}
}
\author{
Ethan Bass
}
