test_that("there are no duplicate peptides after processing", {
  expect_true({
      dat = mspms::prepare_peaks("..//testdata/protein-peptides-lfq_2.csv",
                                 "../testdata/protein-peptides-id_2.csv")
      dat %>%
      dplyr::pull(Peptide) %>%
      dplyr::n_distinct() == nrow(dat)
    })

})
