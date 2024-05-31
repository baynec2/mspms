test_that("works with standard number of AAs)", {
  expect_no_error(mspms::calc_AA_count_of_motif("XXXXXXXX"))
})



test_that("works with non standard number of AAs (so long as even number)", {
  expect_no_error(mspms::calc_AA_count_of_motif("XXXXXXXXXX"))
})


test_that("does not work with odd number of AAs", {
  expect_error(mspms::calc_AA_count_of_motif("XXX"))
})
