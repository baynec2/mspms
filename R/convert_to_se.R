prepared_data = mspms::peaks_prepared_data
design_matrix = mspms::design_matrix

library(dplyr)

convert_to_qf = function(prepared_data,
                         design_matrix){

  # Which columns contain the quantification data?

  quantCols = 4:ncol(prepared_data)
  
  # Reading as a QFeatures object
  qf = QFeatures::readQFeatures(prepared_data,
                                ecol = quantCols,
                                name = "peptides",
                                fnames = "Peptide")
  
  # Adding the col data
  qf$group = design_matrix$group
  qf$condition = design_matrix$condition
  qf$time = design_matrix$time
  
return(qf)
  
}

DEP::plot_missval(qf)

  
QFeatures::impute(qf,method = "MinProb")
?QFeatures::impute
imputeMethods()
?MsCoreUtils::imputeMethods(
)
?QFeatures::impute
