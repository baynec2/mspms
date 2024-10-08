% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics_helper_functions.R
\name{mspms_log2fc}
\alias{mspms_log2fc}
\title{mspms_log2fc}
\usage{
mspms_log2fc(processed_qf, reference_variable = "time", reference_value = 0)
}
\arguments{
\item{processed_qf}{a QFeatures object with a SummarizedExperiment named
"peptides_norm".}

\item{reference_variable}{the variable to used as a reference (denominator
of log2 fold change).}

\item{reference_value}{the value of the reference variable to use as
the reference}
}
\value{
a tibble with the t test statistics for each peptide within each
group with the supplied value at the supplied variable as reference.
}
\description{
calculates the log2fc for each time point within each condition relative to
a specified value for a specified reference variable.
}
\keyword{internal}
