#' calc_sig_zscores
#' determine which zscores are significant at the given alpha for a matrix of zscores
#' @param zscores = a data frame of zscores
#' @param pval = p value threshold for significance. Default is 0.05
#'
#' @return a tibble of significant zscores
#' @export
#' @examples
#' zscores = data.frame(p1 = c(0,1,2,3,4,5),p2 = c(0,-1,-2,-3,-4,-5))
#' rownames(zscores) = c("A","B","C","D","E","F")#'
#' calc_sig_zscores(zscores)

calc_sig_zscores = function(zscores, pval = 0.05){

  #converting p value to zscore threshold. Divide by two since it is two tailed.
  threshold = qnorm(p = pval/2,lower.tail = FALSE)


  sig_zscores = zscores %>%
    tibble::rownames_to_column("AA") %>%
    tidyr::pivot_longer(2:length(.), names_to = "position", values_to = "zscore") %>%
    dplyr::filter(abs(zscore) > threshold) %>%
    dplyr::mutate(aa_position = paste0(.data$AA, ".", .data$position))

  return(sig_zscores)

}
