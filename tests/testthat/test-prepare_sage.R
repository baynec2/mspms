test_that("prepare_sage works", {
  sage_lfq_filepath <- system.file(
    "extdata/sage_lfq.tsv",
    package = "mspms"
  )
  colData_filepath <- system.file("extdata/colData.csv",
    package = "mspms"
  )

  expect_no_error(prepare_sage(sage_lfq_filepath, colData_filepath))
})
