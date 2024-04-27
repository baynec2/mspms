test_that("there are no duplicate peptides after processing", {
  expect_true({
    dat = mspms::prepare_pd("../testdata/proteome_discoverer_output.xlsx")

    dat %>%
      dplyr::pull(Peptide) %>%
      dplyr::n_distinct() == nrow(dat)
  })

})
