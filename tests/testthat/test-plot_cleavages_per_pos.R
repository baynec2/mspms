test_that("plot_cleavages_per_pos() works when given expected data", {
  sig <- log2fc_t_test_data %>%
    dplyr::filter(p.adj <= 0.05, log2fc > 3)

  p1 <- plot_cleavages_per_pos(sig)
  expect_no_error(p1)
})

test_that("plot_cleavages_per_pos() gives error with unexpected data", {
  sig <- log2fc_t_test_data %>%
    dplyr::filter(p.adj <= 0.05, log2fc > 3) %>%
    dplyr::select(-cleavage_pos)

  expect_error(plot_cleavages_per_pos(sig))
})
