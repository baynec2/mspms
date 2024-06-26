% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_all_cleavages.R
\name{calculate_all_cleavages}
\alias{calculate_all_cleavages}
\title{calculate_all_cleavages}
\usage{
calculate_all_cleavages(peptide_library_seqs, n_AA_after_cleavage = 4)
}
\arguments{
\item{peptide_library_seqs}{= these are the sequences from the peptide
library. Should be the same length}

\item{n_AA_after_cleavage}{= this is the number of AA after cleavage sites
you want to produce peptides for.}
}
\value{
a vector of all the possible cleavages for the peptide
 library sequences
}
\description{
calculate all possible cleavages for a defined peptide library containing
 peptides of the same length.
}
\examples{
calculate_all_cleavages(mspms::peptide_library$library_real_sequence,
  n_AA_after_cleavage = 4
)
}
