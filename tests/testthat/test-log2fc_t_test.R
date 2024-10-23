test_that("log2fc_t_test() works using time as reference ", {
    expect_no_error(log2fc_t_test(mspms::processed_qf,
        "time",
        reference_value = 0
    ))
})
# 
# test_that("log2fc_t_test() thows an error when a value is not in data ", {
#     expect_error(log2fc_t_test(mspms::processed_qf,
#         "condition",
#         reference_value = "not present in data"
#     ))
# })
