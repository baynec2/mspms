#' extract_re
#' 
#' extract regular expressions from an icelogo matrix. 
#'
#' @param icelogo_matrix = this is an icelogo matrix, as produced by 
#' prepare_icelogo_data()
#' @param type = "positive" or "negative" depending on whether you want to look
#' for positively or negatively enriched amino acids
#' @param threshold = the threshold to consider amino acids for inclusion. If 
#' type is positive, will return AAs with values greater than the threshold.
#' If type is negative, will return AAs with values less than the threshold.
#' 
#' If no amino acids are found greater than the threshold, a . will be returned 
#' in regular expression. 
#'
#' @return a character vector containing a regular expression extracted from 
#' the icelogo matrix based on user defined parameters. 
#' @export
#'
#' @examples
#'icelogo_matrix = data.frame(P4 = c(10,NA,NA),
#'                            P3 = c(NA,10,NA),
#'                           P2 = c(NA,NA,10),
#'                           P1 = c(NA,NA,NA),
#'                            `P1'` = c(NA,NA,NA),
#'                            `P2'` = c(10,10,10),
#'                            `P3'` = c(10,10,NA),
#'                            `P4'` = c(NA,10,10))
#'rownames(icelogo_matrix) = c("A","B","C")
#'icelogo_matrix = as.matrix(icelogo_matrix)
#'extract_re(icelogo_matrix)
#' 
extract_re <- function(icelogo_matrix,
                       type = "positive",
                       threshold = 0) {
  ## Defining function to finding row names from icelogo matrix for
  # enriched AAs depending on type and percent threshold cutoff.
  find_rownames <- function(x) {
    if (type == "positive") {
      row_name <- names(which(x > threshold))
    } else if (type == "negative") {
      row_name <- names(which(x < threshold))
    } else {
      stop("type must be either positive or negative")
    }
    if (length(row_name) == 0) {
      row_name <- "."
    }
    return(row_name)
  }
  ## Defining function to extract regular expressions from row names
  re <- function(row_names) {
    ## Preparing regular expressions from Row names)
    step1 <- purrr::map(row_names, paste, collapse = "|")
    step2 <- purrr::map(step1, append, "(", after = 0)
    step3 <- purrr::map(step2, append, ")")
    step4 <- purrr::map(step3, paste, collapse = "")
    step5 <- unlist(step4)
    step6 <- paste(step5, collapse = "")
    return(step6)
  }

  rn <- apply(icelogo_matrix, 2, find_rownames)
  regular_expression <- re(rn)

  return(regular_expression)
}
