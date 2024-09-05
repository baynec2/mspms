test_that("t-test gives expected means", {

    summarized = mspms::mspms_data %>% 
      group_by(Peptide, Condition, Time) %>% 
      summarise(mean =mean(value))
  
  t = mspms::mspms_data
})
summarized = mspms::mspms_data %>% 
  dplyr::group_by(Peptide, condition, time) %>% 
  dplyr::summarise(mean = mean(value))

t = mspms::mspms_t_tests(mspms::mspms_data)
