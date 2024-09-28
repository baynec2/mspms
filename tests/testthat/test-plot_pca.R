test_that("plot_pca does not generate an error when testing with expected data", {
    expect_no_error(plot_pca(mspms::mspms_tidy_data))
})
