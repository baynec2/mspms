test_that("plot_volcano() does not generate an error with expected data", {
  expect_no_error(plot_volcano(mspms::log2fc_t_test_data))
})
