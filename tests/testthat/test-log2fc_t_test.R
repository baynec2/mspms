test_that("log2fc_t_test gives expected means", {
  # Had an instance where the means didn't match what was expected
  # Testing here. 
  mspms_data <- mspms::mspms_data %>% 
    mutate(time = as.character(time))
  
  stats <- mspms::log2fc_t_test(mspms_data)
  
  expected_control <- mspms_data %>% 
    dplyr::filter(time == 0) %>% 
    dplyr::group_by(Peptide,condition) %>% 
    dplyr::summarise(expected_control_mean = mean(value,na.rm = TRUE))
  
  expected_sample <- mspms_data %>% 
    dplyr::group_by(Peptide,condition,time) %>% 
    dplyr::summarise(expected_sample_mean = mean(value,na.rm = TRUE))
  
  stats_control = dplyr::left_join(stats,
                                  expected_control,
                                  by = c("Peptide","condition"))
  
  all = dplyr::left_join(stats_control,
                         expected_sample,
                         by = c("Peptide","condition","time"))
  expect_equal(all$expected_control_mean, all$control_mean)
  expect_equal(all$expected_sample_mean, all$sample_mean)
})
