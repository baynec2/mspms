test_that("there are no duplicate peptides after processing", {
  expect_true({
    dat = mspms::prepare_pd("../testdata/proteome_discovere_output.xlsx")

    dat %>%
      dplyr::pull(Peptide) %>%
      dplyr::n_distinct() == nrow(dat)
  })

})
