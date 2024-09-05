test_that("log2fc works as expected", {
  # Had an instance where the means didn't match what was expected
  # Testing here. 
  log2fc <- mspms::mspms_log2fc(mspms::mspms_data)
  
  expected_control <- mspms::mspms_data %>% 
    dplyr::filter(time == 0) %>% 
    dplyr::group_by(Peptide,condition) %>% 
    dplyr::summarise(expected_control_mean = mean(value,na.rm = TRUE))
  
  expected_sample <- mspms::mspms_data %>% 
    dplyr::group_by(Peptide,condition,time) %>% 
    dplyr::summarise(expected_sample_mean = mean(value,na.rm = TRUE))
  
  l2fc_control = dplyr::left_join(log2fc,
                         expected_control,
                         by = c("Peptide","condition"))
  
  all = dplyr::left_join(l2fc_control,
                                  expected_sample,
                                  by = c("Peptide","condition","time"))
  expect_equal(all$expected_control_mean, all$control_mean)
  expect_equal(all$expected_sample_mean, all$sample_mean)
})


