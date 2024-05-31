test_that("throws an error without row names", {
  expect_error(mspms::calc_AA_fc(matrix(c(0.4,0.2),nrow = 1),
                                 matrix(c(0.2,0.4),nrow = 1)))
})
