% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rmph}
\alias{rmph}
\title{Samples \code{n} outcomes of the multivariate phase-type distribution}
\usage{
rmph(n, phase_type_graph, rewards)
}
\arguments{
\item{n}{The number of random samples to draw}

\item{phase_type_graph}{A reference to the graph created by \code{\link[=create_graph]{create_graph()}}}

\item{rewards}{rewards which should be applied to the phase-type distribution. Must have rows equal to \code{\link[=vertices_length]{vertices_length()}} and columns the dimensions of the multivariate phase-type distribution}
}
\value{
A matrix of numeric values, with the parameter \code{rewards} number of rows
}
\description{
Samples \code{n} outcomes of the multivariate phase-type distribution
}
\seealso{
\code{\link[=rph]{rph()}}
}
