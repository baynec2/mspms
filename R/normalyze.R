#' normalyze: this function normalyzes the data using the NormalyzerDE package. 
#'
#' @param prepared_data = this is the data that has been prepared
#' @param design_matrix = this is the design matrix
#' @param outdir = this is the directory you would like NormalyzerDE::nomralyzer to write output to. 
#'
#' @return
#' files written to the outdir specified
#' the normalyzed data is then read in as a data frame. 
#' @export
#'
#' @examples
normalyze = function(prepared_data,design_matrix,outdir = getwd()){
  
  # Extracting only the data
  dataOnly = prepared_data[, design_matrix$sample]

  #Extracting only the annotation columns
  annotOnly = prepared_data[, !(colnames(prepared_data) %in% design_matrix$sample)]
  
  # Making a summarized Experiment object
  sumExpObj = SummarizedExperiment::SummarizedExperiment(
    as.matrix(dataOnly),
    colData = design_matrix,
    rowData = annotOnly)
  

  # Making the jobName be the date folowed by brief tag
  jobName = paste0(Sys.Date(),"_msp-ms_normalyze_output")
  ## Performing the normalization, this function outputs a bunch of stuff to a directory, including plots
  NormalyzerDE::normalyzer(jobName= jobName,
                           experimentObj = sumExpObj,
                           outdir)
  
  ## Reading in the data frame we just produced that we care about ##
  norm_data = readr::read_delim(paste0(outdir,"/",jobName,"/","median-normalized.txt")) %>% 
    #Performing a reverse log2 transformation
    dplyr::mutate(across(30:length(.),~2**.x))

  return(norm_data)
}

