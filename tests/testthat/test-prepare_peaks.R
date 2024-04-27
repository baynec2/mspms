test_that("there are no duplicate peptides after processing", {
  expect_true({
      dat = mspms::prepare_peaks(system.file("tests/testdata/protein-peptides-lfq_2.csv",package = "mspms"),
                                 system.file("tests/testdata/protein-peptides-id_2.csv",package = "mspms"))
      dat %>%
      dplyr::pull(Peptide) %>%
      dplyr::n_distinct() == nrow(dat)
    })

})
