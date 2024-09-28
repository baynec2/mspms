test_that("plot_nd_peptides() does not generate an error when data is as
          expected", {
    expect_no_error(plot_nd_peptides(mspms::processed_qf))
})
