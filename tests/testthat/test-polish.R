test_that("all peptide cleavages are dropped", {
  
  cleavage_added_data = tibble::tibble(nterm = "AAAAAAAA",
         cterm = "AAAAAAAA",
         cterm_cleavage_pos = 13,
         nterm_cleavage_pos = 1) 
  
  polished = mspms::polish(cleavage_added_data)
  
  expect_equal(nrow(polished),0)
})

test_that("cterm cleavages retained", {
  
  cleavage_added_data = tibble::tibble(nterm = NA,
                                       cterm = "AAAAAAAA",
                                       cterm_cleavage_pos = 13,
                                       nterm_cleavage_pos = NA) 
  
  polished = mspms::polish(cleavage_added_data)
  
  expect_equal(nrow(polished),1)
})

test_that("nterm retained", {
  
  cleavage_added_data = tibble::tibble(nterm = "AAAAAAAA",
                                       cterm = NA,
                                       cterm_cleavage_pos = NA,
                                       nterm_cleavage_pos = 1) 
  
  polished = mspms::polish(cleavage_added_data)
  
  expect_equal(nrow(polished),1)
})

test_that("full length peptides retained", {
  
  cleavage_added_data = tibble::tibble(nterm = NA,
                                       cterm = NA,
                                       cterm_cleavage_pos = NA,
                                       nterm_cleavage_pos = NA) 
  
  polished = mspms::polish(cleavage_added_data)
  
  expect_equal(nrow(polished),1)
})

