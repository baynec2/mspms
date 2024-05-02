test_that("there are no duplicate peptides after processing", {
  expect_true({
      dat = mspms::prepare_peaks("../testdata/protein-peptides-lfq_2.csv",
                                 "../testdata/protein-peptides-id_2.csv")
      dat %>%
      dplyr::pull(Peptide) %>%
      dplyr::n_distinct() == nrow(dat)
    })

})

test_that("peaks lfq column error works", {
    expect_error({
    mspms::prepare_peaks("../testdata/protein-peptides-lfq_bad_names.csv",
                         "../testdata/protein-peptides-id.csv")

  })

})


test_that("peaks id column error works", {
  expect_error({
    mspms::prepare_peaks("../testdata/protein-peptides-lfq.csv",
                          "../testdata/protein-peptides-id_bad_names.csv")
  })

})

