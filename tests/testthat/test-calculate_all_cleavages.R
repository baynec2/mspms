test_that("calculate all cleavages works", {
  expect_no_error(calculate_all_cleavages(mspms::peptide_library$library_real_sequence))
})
