test_that("limma_stats works with one condition", {
  expect_no_error(mspms::limma_stats(mspms::processed_qf))
})

