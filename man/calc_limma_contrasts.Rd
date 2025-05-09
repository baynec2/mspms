% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics_helper_functions.R
\name{calc_limma_contrasts}
\alias{calc_limma_contrasts}
\title{calc_limma_contrasts}
\usage{
calc_limma_contrasts(colData, design_mat)
}
\arguments{
\item{colData}{colData from mspms experiment}

\item{design_mat}{design_mat as returned by calc_limma_design_matrix}
}
\value{
a contrast matrix
}
\description{
Calculates limma contrasts given colData. The contrasts returned are pairwise
relative to T0 for each timepoint assayed.
}
\keyword{internal}
