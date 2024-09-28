# Testing for case where there is no cleavage site
testthat::test_that("works with no cleavage site", {
    testthat::expect_equal(
        mspms:::cterm_cleavage(
            "ABCDEFGHIJKLMN", "ABCDEFGHIJKLMN",
            "ABCDEFGHIJKLMN", 4
        ),
        tibble::tibble(
            peptide = "ABCDEFGHIJKLMN",
            cterm = NA,
            cterm_cleavage_pos = NA
        )
    )
})


# Testing for case where there is no cleavage site
testthat::test_that("works with one clevage site", {
    testthat::expect_equal(
        mspms:::cterm_cleavage(
            "ABCDEFGHIJKLM_N", "ABCDEFGHIJKLMN",
            "ABCDEFGHIJKLMN", 4
        ),
        tibble::tibble(
            peptide = "ABCDEFGHIJKLM_N",
            cterm = "JKLMNXXX",
            cterm_cleavage_pos = 13
        )
    )
})


# Testing for case where there are two cleavage sites
testthat::test_that("works with two clevage site", {
    testthat::expect_equal(
        mspms:::cterm_cleavage(
            "A_BCDEFGHIJKLM_N", "ABCDEFGHIJKLMN",
            "ABCDEFGHIJKLMN", 4
        ),
        tibble::tibble(
            peptide = "A_BCDEFGHIJKLM_N",
            cterm = "JKLMNXXX",
            cterm_cleavage_pos = 13
        )
    )
})


# Testing for case where we look for 6 AAs past the cleavage site
testthat::test_that("works with 6 AA past", {
    testthat::expect_equal(
        mspms:::cterm_cleavage(
            "A_BCDEFGHIJKLM_N", "ABCDEFGHIJKLMN",
            "ABCDEFGHIJKLMN", 6
        ),
        tibble::tibble(
            peptide = "A_BCDEFGHIJKLM_N",
            cterm = "HIJKLMNXXXXX",
            cterm_cleavage_pos = 13
        )
    )
})


# Testing for extreme case
testthat::test_that("works with 10 AA past", {
    testthat::expect_equal(
        mspms:::cterm_cleavage(
            "A_BCDEFGHIJKLM_N", "ABCDEFGHIJKLMN",
            "ABCDEFGHIJKLMN", 10
        ),
        tibble::tibble(
            peptide = "A_BCDEFGHIJKLM_N",
            cterm = "DEFGHIJKLMNXXXXXXXXX",
            cterm_cleavage_pos = 13
        )
    )
})
