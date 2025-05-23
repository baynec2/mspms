% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_diann.R
\name{prepare_diann}
\alias{prepare_diann}
\title{prepare_diann}
\usage{
prepare_diann(
  precursor_filepath,
  colData_filepath,
  peptide_library = mspms::peptide_library,
  n_residues = 4
)
}
\arguments{
\item{precursor_filepath}{filepath to report.pr_matrix.tsv file exported
from DIA-NN.}

\item{colData_filepath}{file path to .csv file containing colData.
Must have columns named "quantCols","group","condition",and "time".}

\item{peptide_library}{peptide library used with experiment. Contains
columns "library_id", "library_match_sequence", and "library_real_sequence".}

\item{n_residues}{the number of amino acid residues before and after the
cleavage site to generate a cleavage seq for.}
}
\value{
a QFeatures object.
}
\description{
prepare data from the pr_matrix.tsv diann output. This can be either from
DIA-NN or from Fragpipe (as it uses DIA-NN for quantification internally for
MSFragger-DIA workflows)
}
\examples{
precursor_filepath <- system.file(
  "extdata/diann_report.pr_matrix.tsv",
  package = "mspms"
)
colData_filepath <- system.file("extdata/diann_colData.csv",package = "mspms")
prepare_diann(precursor_filepath,colData_filepath)
}
