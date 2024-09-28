test_that("calculate all cleavages works", {
    expect_no_error(
        mspms::calculate_all_cleavages(
            mspms::peptide_library$library_real_sequence
        )
    )
})

test_that("calculate all cleavages generates cleavage motifs of the same
          length", {
    background_universe <- mspms::calculate_all_cleavages(
        mspms::peptide_library$library_real_sequence
    )
    nchar <- unique(sapply(background_universe, nchar))
    expect_equal(
        length(nchar),
        1
    )
})
