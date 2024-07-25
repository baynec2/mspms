test_that("two plots are in the file", {
  
  qc_check_data = tibble::tibble(perc_undigested_library_detected = runif(4,0,100),
                           perc_library_detected = runif(4,0,100),
                           group = rep(c("A","B"),2))
  
  p = plot_qc_check(qc_check_data)
  
  expect_length(p,2)

})
