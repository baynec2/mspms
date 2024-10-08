% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_helper_functions.R
\name{prepare_for_PCA}
\alias{prepare_for_PCA}
\title{prepare_for_PCA()}
\usage{
prepare_for_PCA(mspms_tidy_data, value_colname = "peptides_norm")
}
\arguments{
\item{mspms_tidy_data}{tidy mspms data (prepared from QFeatures object
by mspms_tidy())}

\item{value_colname}{the name of the column containing values.}
}
\value{
a tibble
}
\description{
prepare QFeatures object for PCA analysis
}
\keyword{internal}
