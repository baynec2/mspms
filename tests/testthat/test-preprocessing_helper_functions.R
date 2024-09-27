### check_file_is_valid_peaks()###
test_that("check_file_is_valid_peaks() errors when PTM is provided.", {
  peaks_ptm <- readr::read_csv(system.file(
    "extdata/peaks_protein-peptides-lfq.csv",
    package = "mspms"
  )) %>%
    dplyr::mutate(PTM = "PTM")
  expect_error(mspms:::check_file_is_valid_peaks(peaks_ptm))
})

### check_file_is_valid_fragpipe()###
test_that("check_file_is_valid_fragpipe() errors when peaks file is
          provided.", {
  peaks_ptm <- readr::read_csv("tests/testdata/
                               check_file_is_valid_peaks_ptm.csv")
  expect_error(mspms:::check_file_is_valid_fragpipe(peaks_ptm))
})

### consolidate_cleavages()###

test_that("consolidate_cleavages() removes double cleaved peptides", {
  test <- tibble::tibble(
    peptide = "P_AAAAAAAAAAAAA",
    nterm = "P_AAAAAAAA_P",
    nterm_cleavage_pos = 5,
    cterm = "P_AAAAAAA_P",
    cterm_cleavage_pos = 5
  )
  nrow <- nrow(mspms:::consolidate_cleavages(test))
  expect_equal(nrow, 0)
})

### prepared_to_qf ###
test_that("prepared_to_qf() generates an error when supplied a QFeatures
          object", {
  expect_error(mspms:::prepared_to_qf(mspms::peaks_prepared_data,
    colData = mspms::colData
  ))
})


test_that("prepared_to_qf() fails with unexpected colData", {
  
  prepared_data = mspms::processed_qf %>% 
    mspms_tidy("peptides") %>% 
    dplyr::select(peptide,library_id,quantCols,peptides) %>% 
    tidyr::pivot_wider(names_from = "quantCols", values_from = "peptides")
 
   wrong_col_data = mspms::colData %>% 
    mutate(quantCols = 1:length(quantCols))
   
   expect_error(mspms:::prepared_to_qf(prepared_data,wrong_col_data))
    
})

test_that("prepared_to_qf() works with expected colData", {
  
  prepared_data = mspms::processed_qf %>% 
    mspms_tidy("peptides") %>% 
    dplyr::select(peptide,library_id,quantCols,peptides) %>% 
    tidyr::pivot_wider(names_from = "quantCols", values_from = "peptides")

  expect_no_error(mspms:::prepared_to_qf(prepared_data,mspms::colData))
  
})

### load_colData() ###
test_that("load_colData() generates expected expected error when colData is
          off", {
  expect_error(
    mspms:::load_colData(
      "tests/testdata/missing_time.csv",
      "colData must have columns named \"quantCols\",\"group\",
                       \"condition\", and \"time\" you are missing time"
    )
  )
})
