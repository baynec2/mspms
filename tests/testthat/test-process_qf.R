test_that("process_qf() produces the expected assays", {
  processed_qf <- process_qf(mspms::peaks_prepared_data)
  expected_assays <- c(
    "peptides", "peptides_log", "peptides_log_norm",
    "peptides_log_impute_norm", "peptides_norm"
  )
  sum <- sum(names(processed_qf) %!in% expected_assays)
  expect_equal(sum, 0)
})
