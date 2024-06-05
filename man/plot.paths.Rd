% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paths.R
\name{plot.paths}
\alias{plot.paths}
\title{Visualize the paths of a SDTree or SDForest}
\usage{
\method{plot}{paths}(object, plotly = F, selection = NULL, log_scale = F)
}
\arguments{
\item{object}{A paths object}

\item{plotly}{If TRUE the plot is interactive using plotly. Might be slow for large data.}

\item{selection}{A vector of indices of the covariates to be plotted. 
Can be used to plot only a subset of the covariates in case of many covariates.}

\item{log_scale}{If TRUE the y-axis is on a log scale.}
}
\value{
A ggplot object
}
\description{
This function visualizes the variable importance of a SDTree or SDForest
for different complexity parameters. Both the regularization path and
the stability selection path can be visualized.
}
\seealso{
\code{\link{regPath}} \code{\link{stabilitySelection}}
}
\author{
Markus Ulmer
}