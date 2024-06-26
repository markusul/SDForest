% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partDependence.R
\name{plot.partDependence}
\alias{plot.partDependence}
\title{Plot partial dependence}
\usage{
\method{plot}{partDependence}(x, n_examples = 19, ...)
}
\arguments{
\item{x}{An object of class \code{partDependence} returned by \code{\link{partDependence}}.}

\item{n_examples}{Number of examples to plot in addition to the average prediction.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
A ggplot object.
}
\description{
This function plots the partial dependence of a model on a single variable.
}
\seealso{
\code{\link{partDependence}}
}
\author{
Markus Ulmer
}
