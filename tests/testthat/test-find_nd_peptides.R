test_that("function works", {
  expect_no_error(find_nd_peptides(mspms::peaks_prepared_data,
                                   mspms::peptide_library,
                                   mspms::design_matrix))
})
