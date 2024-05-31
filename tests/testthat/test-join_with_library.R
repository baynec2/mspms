test_that("join with library  works", {
  expect_no_error(join_with_library(mspms::imputed, mspms::peptide_library))
})
