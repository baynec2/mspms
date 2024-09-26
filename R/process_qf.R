#' process_qf
#'
#' @param prepared_qf this is a QFeatures object containing a
#' SummarizedExperiment named "peptides"
#'
#' @return a QFeatures object containing a
#' SummarizedExperiments named "peptides","peptides_log","peptides_log_norm",
#' "peptides_log_impute_norm",and "peptides_norm"
#' @export
#' @examples
#' processed_qf <- process_qf(mspms::peaks_prepared_data)
process_qf <- function(prepared_qf) {
  # Log2transform the data, add as assay
  prepared_qf_norm <- QFeatures::addAssay(prepared_qf,
    QFeatures::logTransform(
      prepared_qf[["peptides"]]
    ),
    name = "peptides_log"
  )
  # Normalize the data, add as assay
  prepared_qf_norm <- QFeatures::normalize(prepared_qf_norm,
    i = "peptides_log",
    name = "peptides_log_norm",
    method = "center.median"
  )
  # Impute the data, add as assay
  QF_mod <- QFeatures::impute(prepared_qf_norm,
    method = "QRILC",
    i = "peptides_log_norm",
    name = "peptides_log_impute_norm", MARGIN = 2
  )
  # reverse log2 transform the data
  QF_mod <- QFeatures::aggregateFeatures(QF_mod,
    i = "peptides_log_impute_norm",
    "peptide", name = "peptides_norm",
    fun = rlog2
  )
  return(QF_mod)
}
