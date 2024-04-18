# Testing for case where there is no cleavage site
testthat::test_that("works with no cleavage site", {
  testthat::expect_equal(
    mspms::nterm_cleavage("ABCDEFGHIJKLMN", "ABCDEFGHIJKLMN", "ABCDEFGHIJKLMN",4),
    data.frame(peptide = "ABCDEFGHIJKLMN",
               nterm = NA,
               nterm_cleavage_pos = NA)
  )
})


# Testing for case where there is no cleavage site
testthat::test_that("works with one clevage site", {
  testthat::expect_equal(
    mspms::nterm_cleavage("A_BCDEFGHIJKLM_N", "ABCDEFGHIJKLMN", "ABCDEFGHIJKLMN",4),
    data.frame(peptide = "A_BCDEFGHIJKLM_N",
               nterm = "XXXABCDE",
               nterm_cleavage_pos = 1)
  )
})


# Testing for case where there are two cleavage sites
testthat::test_that("works with two clevage site", {
  testthat::expect_equal(
    mspms::nterm_cleavage("A_BCDEFGHIJKLM_N", "ABCDEFGHIJKLMN", "ABCDEFGHIJKLMN",4),
    data.frame(peptide = "A_BCDEFGHIJKLM_N",
               nterm = "XXXABCDE",
               nterm_cleavage_pos = 1)
  )
})


# Testing for case where we look for 6 AAs past the cleavage site
testthat::test_that("works with 6 AA past", {
  testthat::expect_equal(
    mspms::nterm_cleavage("A_BCDEFGHIJKLM_N", "ABCDEFGHIJKLMN", "ABCDEFGHIJKLMN",6),
    data.frame(peptide = "A_BCDEFGHIJKLM_N",
               nterm = "XXXXXABCDEFG",
               nterm_cleavage_pos = 1)
  )
})


# Testing for extreme case
testthat::test_that("works with 10 AA past", {
  testthat::expect_equal(
    mspms::nterm_cleavage("A_BCDEFGHIJKLM_N", "ABCDEFGHIJKLMN", "ABCDEFGHIJKLMN",10),
    data.frame(peptide = "A_BCDEFGHIJKLM_N",
               nterm = "XXXXXXXXXABCDEFGHIJK",
               nterm_cleavage_pos = 1)
  )
})
