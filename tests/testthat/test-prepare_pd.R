test_that("there are no duplicate peptides after processing", {
  expect_true({
    dat = mspms::prepare_pd(system.file("tests/testdata/proteome_discoverer_output.xlsx",package = "mspms"))

    dat %>%
      dplyr::pull(Peptide) %>%
      dplyr::n_distinct() == nrow(dat)
  })

})
