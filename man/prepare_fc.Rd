% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iceLogo_helper_fuctions.R
\name{prepare_fc}
\alias{prepare_fc}
\title{prepare_fc}
\usage{
prepare_fc(fold_change, sig_zscores)
}
\arguments{
\item{fold_change}{a matrix of the fold changes of the AA by
position.}

\item{sig_zscores}{a tibble of the significant zscores.}
}
\value{
a matrix of the fold changes of the significant AAs at each position.
}
\description{
Prepare fold changes of amino acids by position for Icelogo visualization.
}
\keyword{internal}
