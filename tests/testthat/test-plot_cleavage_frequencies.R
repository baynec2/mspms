test_that("plot_cleavage_frequencies works", {
  expect_no_error(plot_cleavage_frequencies(rep("AAAAAAAA",20)))
})
test_that("plot_cleavage_frequencies doesn't work with mismatched lengths", {
  expect_error(plot_cleavage_frequencies(rep("AAAAAAA",20)))
})
