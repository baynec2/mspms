test_that("prepare_diann works", {
  expect_no_error({
    precursor_filepath <- system.file("extdata/diann_report.pr_matrix.tsv",
                                      package = "mspms")
    colData_filepath <- system.file("extdata/diann_colData.csv",package = "mspms")
    prepare_diann(precursor_filepath,colData_filepath)
  })
})
