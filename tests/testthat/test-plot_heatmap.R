test_that("plot_heatmap() produces no error with expected data", {
  expect_no_error(plot_heatmap(mspms::mspms_tidy_data))
})

test_that("plot_heatmap() produces error with bad argument", {
  expect_error(plot_heatmap(mspms::mspms_tidy_data, scale = "wrong_arg"))
})
