test_that("plot_all_icelogos() produces ggarrange s3 class", {
    sig <- mspms::log2fc_t_test_data %>%
        dplyr::filter(p.adj <= 0.05, log2fc > 3)
    t <- plot_all_icelogos(sig)
    expect_s3_class(t, "ggarrange")
})

test_that("plot_all_icelogos() works with nresidues = 4", {
  sig <- mspms::log2fc_t_test_data %>%
    dplyr::filter(p.adj <= 0.05, log2fc > 3)
    expect_no_error(plot_all_icelogos(sig))
})


test_that("plot_all_icelogos() works with nresidues = 8", {
    
  all_16mers = calculate_all_cleavages(
    mspms::peptide_library$library_real_sequence,n_AA_after_cleavage = 8
    )
  
  peaks_prepared_data = mspms::prepare_peaks(
    system.file("extdata/peaks_protein-peptides-lfq.csv",
                package = "mspms"
    ),
    system.file("extdata/colData.csv", package = "mspms"),
    n_residues = 8
  ) %>% 
    process_qf() %>% 
    log2fc_t_test()
  
   sig <- peaks_prepared_data %>%
    dplyr::filter(p.adj <= 0.05, log2fc > 3)
   
  expect_no_error(plot_all_icelogos(sig,background_universe = all_16mers))
})
