#' rlog2
#' Reverse log2 transformation
#' @param x a numeric value
#'
#' @return a reverse log2 transformed value
#' @keywords internal
rlog2 <- function(x) {
  rlog2 <- 2^x
  return(rlog2)
}
