test_that("there are no duplicate peptides after processing", {
  expect_true({
    dat = mspms::prepare_pd("../testdata/proteome_discoverer_output.xlsx")

    dat %>%
      dplyr::pull(Peptide) %>%
      dplyr::n_distinct() == nrow(dat)
  })

})

test_that("error to find unexpected id columns works", {
  expect_error({
    mspms::prepare_pd(filepath = "../testdata/proteome_discoverer_output_bad_names.xlsx")
  })

})
