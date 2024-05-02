test_that("finds expected cleavage", {
  expect_true(check_for_cleavage_after_AA("GLYWAFAA",c("W","L","Y","F")))
})

test_that("does not find unexpected cleavage", {
  expect_false(check_for_cleavage_after_AA("GLYAAFAA",c("W","L","Y","F")))
})
