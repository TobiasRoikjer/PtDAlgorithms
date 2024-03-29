% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{moments}
\alias{moments}
\title{Computes the first \code{k} moments of the phase-type distribution}
\usage{
moments(phase_type_graph, power, rewards = NULL)
}
\arguments{
\item{phase_type_graph}{A reference to the graph created by \code{\link[=create_graph]{create_graph()}}}

\item{power}{An integer of the first \code{k} moments.}

\item{rewards}{Optional rewards, which should be applied to the phase-type distribution. Must have length equal to \code{\link[=vertices_length]{vertices_length()}}}
}
\value{
A numeric vector of the first \code{k} moments. The first entry is the first moment (mean)
}
\description{
This function invokes \code{\link[=expected_waiting_times]{expected_waiting_times()}} consequtively to find the first moments,
given by the \code{power} argument
}
\examples{
graph <- ptdalgorithms::create_graph(4)
v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 1)
ptdalgorithms::add_edge(v1, v2, 4)
ptdalgorithms::add_edge(v2, a, 10)
ptdalgorithms::moments(graph, 3) # =>
  (0.350000 0.097500 0.025375)
ptdalgorithms::moments(graph, 3, c(0,2,1,0)) # =>
  (0.600 0.160 0.041)
}
\seealso{
\code{\link[=expected_waiting_time]{expected_waiting_time()}}

\code{\link[=expectation]{expectation()}}

\code{\link[=variance]{variance()}}

\code{\link[=covariance]{covariance()}}
}
