test_that("function works", {
  prepared_data = tibble::tibble(Peptide = c("A_AAAAAAA","AAAAAAAAAAAA"),
                         RT = c(1,2),
                         `Protein Accession` = c(1,1),
                         sample_1 = c(100,100))
  
  peptide_library = tibble::tibble(library_reference_id = 1)
  
  design_matrix = tibble::tibble(sample = "sample_1",time = 0, condition = 0)
  
  out = qc_check(prepared_data, peptide_library, design_matrix)
  
  expect_equal(out$perc_library_detected, 100)
  expect_equal(out$perc_undigested_library_detected, 100)
})
