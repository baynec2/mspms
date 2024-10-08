% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{processed_qf}
\alias{processed_qf}
\title{processed_qf
A QFeatures object prepared from PEAKS data of Cathepsin data that has been
processed (imputation/normalization)}
\format{
## `peaks_prepared_data`
An instance of class QFeatures containing 5 assays:
[1] peptides: SummarizedExperiment with 2071 rows and 42 columns
[2] peptides_log: SummarizedExperiment with 2071 rows and 42 columns
[3] peptides_log_norm: SummarizedExperiment with 2071 rows and 42 columns
[4] peptides_log_impute_norm: SummarizedExperiment with 2071 rows and 42
columns
[5] peptides_norm: SummarizedExperiment with 2071 rows and 42 columns
\describe{
  \item{peptides}{Peptide Sequence Detected}
  ...
}
}
\source{
<mspms processed data originally from PEAKS files found in
 "tests/testdata/protein-peptides-id.csv" and
 "tests/testdata/protein-peptides-lfq.csv">
}
\usage{
processed_qf
}
\description{
processed_qf
A QFeatures object prepared from PEAKS data of Cathepsin data that has been
processed (imputation/normalization)
}
\keyword{datasets}
