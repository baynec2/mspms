test_that("prepare_peaks() works as expected.", {
    expect_no_error({
        mspms::prepare_peaks(
            system.file("extdata/peaks_protein-peptides-lfq.csv",
                package = "mspms"
            ),
            system.file("extdata/colData.csv", package = "mspms")
        )
    })
})
