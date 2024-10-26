test_that("plot_qc_check() does not give an error when expected data as
          expected is supplied", {
  expect_no_error(plot_qc_check(mspms::processed_qf))
})
