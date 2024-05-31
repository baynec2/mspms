test_that("imputation works", {
  expect_no_error(impute(mspms::outliers, noise = 0.05, nsd = 1))
})
