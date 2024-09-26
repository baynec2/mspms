test_that("plot_time_course() does not generate an error when data is as
          expected", {
  expect_no_error(plot_time_course(mspms::mspms_tidy_data))
})
