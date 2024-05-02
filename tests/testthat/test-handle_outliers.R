test_that("handle outliers works", {
  expect_no_error(handle_outliers(mspms::normalyzed_data,mspms::design_matrix))

})
