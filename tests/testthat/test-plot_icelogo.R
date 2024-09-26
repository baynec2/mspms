test_that("plot_icelogo() does not generate an error with expected data", {
  cleavage_seqs <- mspms::log2fc_t_test_data %>%
    dplyr::filter(
      p.adj <= 0.05, log2fc > 3,
      condition == "CatA"
    ) %>%
    dplyr::pull(cleavage_seq)

  expect_no_error(plot_icelogo(cleavage_seqs))
})
