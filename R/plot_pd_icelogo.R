#' plot_pd_icelogo
#'
#' plot a icelogo of the percent difference.
#'
#' @param prepared_sig_pdiff_data = this is the prepared significant percent difference matrix.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' prepared_sig_pdiff_data = data.frame(P1 = c(50,10),P2 = c(-50,-10))
#' rownames(prepared_sig_pdiff_data) = c('R','K')
#' prepared_sig_pdiff_data = as.matrix(prepared_sig_pdiff_data)
#' plot_pd_icelogo(prepared_sig_pdiff_data)

plot_pd_icelogo = function(prepared_sig_pdiff_data){

  nchar_cleav = ncol(prepared_sig_pdiff_data)

  p1 = ggseqlogo::ggseqlogo(prepared_sig_pdiff_data,
                            font = "helvetica_light" ,
                            method = 'custom',
                            seq_type = "AA") +
    ggplot2::scale_x_continuous(labels = paste0("P", c((nchar_cleav / 2):1,
                                                       paste0(
                                                         seq_len(nchar_cleav / 2), "'"
                                                       ))),
                                breaks = seq_len(nchar_cleav))+
    ggplot2::ylab("Percent Difference")

  return(p1)
}
