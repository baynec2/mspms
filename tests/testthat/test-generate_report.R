test_that("generate report does not generate an error", {
  expect_no_error(
    generate_report(
      mspms::peaks_prepared_data,
      mspms::peptide_library,
      n = 4
    )
  )
})
