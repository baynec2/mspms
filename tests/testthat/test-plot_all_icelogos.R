test_that("plot_all_icelogos() produces ggarrange s3 class", {
  sig <- mspms::log2fc_t_test_data %>%
    dplyr::filter(p.adj <= 0.05, log2fc > 3)
  t <- plot_all_icelogos(sig)
  expect_s3_class(t, "ggarrange")
})
