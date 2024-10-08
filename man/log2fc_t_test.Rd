% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log2fc_t_test.R
\name{log2fc_t_test}
\alias{log2fc_t_test}
\title{log2fc_t_test}
\usage{
log2fc_t_test(processed_qf, reference_variable = "time", reference_value = 0)
}
\arguments{
\item{processed_qf}{mspms data in a QFeatures object.}

\item{reference_variable}{the colData variable to use as reference}

\item{reference_value}{the value of the colData variable to use as reference}
}
\value{
a tibble containing log2fc and t test statistics
}
\description{
Calculates the log2 fold change and t-test statistics given a user specified
reference variable and value.
}
\examples{
log2fc_and_t_test <- log2fc_t_test(mspms::processed_qf)
}
