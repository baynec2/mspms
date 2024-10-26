### remaining_cd_names() ###
test_that("remaining_cd_names() returns the expected missing name, when time is
provided as reference", {
  time <- mspms:::remaining_cd_names(mspms::processed_qf, "time")
  expect_equal(time, "condition")
})

test_that("remaining_cd_names() returns the expected missing name, when condition
is provided as reference", {
  condition <- mspms:::remaining_cd_names(mspms::processed_qf, "condition")
  expect_equal(condition, "time")
})

test_that("remaining_cd_names() throws error when column not in data is
          attempted to be passed", {
  expect_error(mspms:::remaining_cd_names(mspms::processed_qf, "not_here"))
})

### mspms_log2fc() ###
test_that("mspms_log2fc() returns does not error when data is as expected", {
  expect_no_error(mspms:::mspms_log2fc(
    mspms::processed_qf,
    "time",
    0
  ))
})

test_that("mspms_log2fc() returns does not error when time is supplied as a
          character", {
  expect_no_error(mspms:::mspms_log2fc(
    mspms::processed_qf,
    "time",
    "0"
  ))
})

test_that("mspms_log2fc() errors when a variable not in data is provided", {
  expect_error(mspms:::mspms_log2fc(
    mspms::processed_qf,
    "not_here",
    "0"
  ))
})

### mspms_t_tests() ###
test_that("mspms_t_tests() does not error when data is as expected", {
  expect_no_error(mspms:::mspms_t_tests(
    mspms::processed_qf,
    "time",
    "0"
  ))
})

# test_that("mspms_t_tests() does not error when numeric value of reference
#           variable is provided", {
#     expect_no_error(mspms:::mspms_t_tests(
#         mspms::processed_qf,
#         "time",
#         0
#     ))
# })
#
# test_that("mspms_t_tests() errors when nonexistant value is provided", {
#     expect_error(mspms:::mspms_t_tests(
#         mspms::processed_qf,
#         "not_here",
#         "not_here"
#     ))
# })
