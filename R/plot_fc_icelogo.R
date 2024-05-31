#' plot_fc_icelogo
#'
#' plot an icelogo based on the fold change.
#'
#' @param prepared_fc_data = matrix with prepared fold change data
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' prepared_fc_data = data.frame(P1 = c(5,10),P2 = c(1,-1))
#' rownames(prepared_fc_data) = c('R','K')
#' prepared_fc_data = as.matrix(prepared_fc_data)
#' plot_fc_icelogo(prepared_fc_data)
#'
plot_fc_icelogo = function(prepared_fc_data){
  nchar_cleav = ncol(prepared_fc_data)

  # Plotting the motif
  p1 = ggseqlogo::ggseqlogo(prepared_fc_data,
                            font = "helvetica_light" ,
                            method = 'custom',
                            seq_type = "AA") +
    ggplot2::scale_x_continuous(labels = paste0("P", c((nchar_cleav / 2):1,
                                                       paste0(
                                                         seq_len(nchar_cleav / 2), "'"
                                                       ))),
                                breaks = seq_len(nchar_cleav))+
    ggplot2::ylab("fold change")

  return(p1)
}
