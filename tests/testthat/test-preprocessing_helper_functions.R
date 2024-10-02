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


### Testing nterm cleavages ###

test_that("nterm_cleavages() works as expected with n_residue = 1", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 1)
  
  expect_equal(out$nterm,"GH")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 2", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 2)
  
  expect_equal(out$nterm,"FGHI")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 3", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 3)
  
  expect_equal(out$nterm,"EFGHID")
  expect_equal(out$nterm_cleavage_pos,9)
  
})


test_that("nterm_cleavages() works as expected with n_residue = 4", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
out = mspms:::nterm_cleavage(peptide_sequence,
                       library_match_sequence,
                       library_real_sequence,
                       n_residues = 4)

expect_equal(out$nterm,"YEFGHIDH")
expect_equal(out$nterm_cleavage_pos,9)

})

test_that("nterm_cleavages() works as expected with n_residue = 5", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =  "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 5)
  
  expect_equal(out$nterm,"VYEFGHIDHL")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 6", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =  "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 6)
  
  expect_equal(out$nterm,"TVYEFGHIDHLX")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 7", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =  "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 7)
  
  expect_equal(out$nterm,"ATVYEFGHIDHLXX")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 8", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 8)
  
  expect_equal(out$nterm,"VATVYEFGHIDHLXXX")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 9", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 9)
  
  expect_equal(out$nterm,"LVATVYEFGHIDHLXXXX")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 10", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 10)
  
  expect_equal(out$nterm,"XLVATVYEFGHIDHLXXXXX")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 11", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 11)
  
  expect_equal(out$nterm,"XXLVATVYEFGHIDHLXXXXXX")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 12", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 12)
  
  expect_equal(out$nterm,"XXXLVATVYEFGHIDHLXXXXXXX")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 13", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 13)
  
  expect_equal(out$nterm,"XXXXLVATVYEFGHIDHLXXXXXXXX")
  expect_equal(out$nterm_cleavage_pos,9)
  
})

test_that("nterm_cleavages() works as expected with n_residue = 14", {
  peptide_sequence = "G_HIDHL"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 14)
  
  expect_equal(out$nterm,"XXXXXLVATVYEFGHIDHLXXXXXXXXX")
  expect_equal(out$nterm_cleavage_pos,9)
  
})


test_that("nterm_cleavages() finds correct position", {
  peptide_sequence = "T_VYEFGH"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::nterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 4)
  
  expect_equal(out$nterm,"LVATVYEF")
  expect_equal(out$nterm_cleavage_pos,4)
  
})


### Testing Cterm ###

test_that("cterm_cleavages() works as expected with n_residue = 1", {
  peptide_sequence = "GHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 1)
  
  expect_equal(out$cterm,"HL")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 2", {
  peptide_sequence = "GHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 2)
  
  expect_equal(out$cterm,"DHLX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 3", {
  peptide_sequence = "GHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 3)
  
  expect_equal(out$cterm,"IDHLXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 4", {
  peptide_sequence = "GHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 4)
  
  expect_equal(out$cterm,"HIDHLXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 5", {
  peptide_sequence = "GHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 5)
  
  expect_equal(out$cterm,"GHIDHLXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 6", {
  peptide_sequence = "GHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 6)
  
  expect_equal(out$cterm,"FGHIDHLXXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 7", {
  peptide_sequence = "GHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 7)
  
  expect_equal(out$cterm,"EFGHIDHLXXXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 8", {
  peptide_sequence = "GHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 8)
  
  expect_equal(out$cterm,"YEFGHIDHLXXXXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 9", {
  peptide_sequence = "LVATVYEFGHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 9)
  
  expect_equal(out$cterm,"VYEFGHIDHLXXXXXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})


test_that("cterm_cleavages() works as expected with n_residue = 10", {
  peptide_sequence = "LVATVYEFGHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 10)
  
  expect_equal(out$cterm,"TVYEFGHIDHLXXXXXXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 11", {
  peptide_sequence = "LVATVYEFGHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 11)
  
  expect_equal(out$cterm,"ATVYEFGHIDHLXXXXXXXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 12", {
  peptide_sequence = "LVATVYEFGHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 12)
  
  expect_equal(out$cterm,"VATVYEFGHIDHLXXXXXXXXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 13", {
  peptide_sequence = "LVATVYEFGHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"

  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 13)
  
  expect_equal(out$cterm,"LVATVYEFGHIDHLXXXXXXXXXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() works as expected with n_residue = 14", {
  peptide_sequence = "LVATVYEFGHIDH_L"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"

  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 14)
  
  expect_equal(out$cterm,"XLVATVYEFGHIDHLXXXXXXXXXXXXX")
  expect_equal(out$cterm_cleavage_pos,13)
  
})

test_that("cterm_cleavages() finds correct position", {
  peptide_sequence = "VATV_Y"
  library_match_sequence = "LVATVYEFGHIDHL"
  library_real_sequence =   "LVATVYEFGHIDHL"
  
  out = mspms:::cterm_cleavage(peptide_sequence,
                               library_match_sequence,
                               library_real_sequence,
                               n_residues = 4)
  
  expect_equal(out$cterm,"VATVYEFG")
  expect_equal(out$cterm_cleavage_pos,5)
  
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
    prepared_data <- mspms::processed_qf %>%
        mspms_tidy("peptides") %>%
        dplyr::select(peptide, library_id, quantCols, peptides) %>%
        tidyr::pivot_wider(names_from = "quantCols", values_from = "peptides")

    wrong_col_data <- mspms::colData %>%
        dplyr::mutate(quantCols = 1:length(quantCols))

    expect_error(mspms:::prepared_to_qf(prepared_data, wrong_col_data))
})

test_that("prepared_to_qf() works with expected colData", {
    prepared_data <- mspms::processed_qf %>%
        mspms_tidy("peptides") %>%
        dplyr::select(peptide, library_id, quantCols, peptides) %>%
        tidyr::pivot_wider(names_from = "quantCols", values_from = "peptides")

    expect_no_error(mspms:::prepared_to_qf(prepared_data, mspms::colData))
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
