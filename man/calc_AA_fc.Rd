% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_AA_fc.R
\name{calc_AA_fc}
\alias{calc_AA_fc}
\title{calc_AA_fc}
\usage{
calc_AA_fc(experimental_prop_matrix, background_prop_matrix, sig_zscores)
}
\arguments{
\item{experimental_prop_matrix}{this is a matrix of the experimental
proportions (from your vector of cleavage sequences) at each position.}

\item{background_prop_matrix}{this is a matrix of the background proportions
of AAs at each position}

\item{sig_zscores}{this is a tibble of the significant zscores.}
}
\value{
a matrix
}
\description{
calculate the fold change of each AA by position.
}
\examples{
\dontshow{if (isTRUE(FALSE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontshow{\}) # examplesIf}
}
