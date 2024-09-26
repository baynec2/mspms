### calc_per_samples_library_nd() ###
test_that("all library ids are in output of calc_per_samples_library_nd()", {
  nd_peptides <- mspms:::calc_per_samples_library_nd(
    mspms::peaks_prepared_data,
    peptide_library_ids = mspms::peptide_library$library_id
  )
  li <- unique(nd_peptides$library_id)
  total <- mspms::peptide_library$library_id
  sum <- sum(total %!in% li)
  expect_equal(sum, 0)
})

### prepare_qc_check() ###
test_that("prepare_qc_check() generates percentages 100 and below", {
  out <- mspms:::prepare_qc_check_data(mspms::processed_qf)
  per_library_id_detected <- max(out$per_library_id_detected)
  expect_lt(per_library_id_detected, 100)
})


test_that("prepare_qc_check() does not generate a percentage below 100.", {
  out <- mspms:::prepare_qc_check_data(mspms::processed_qf)
  per_library_id_detected <- min(out$per_library_id_detected)
  expect_gt(per_library_id_detected, 0)
})

### icelogo_col_scheme() ###
test_that("icelogo_col_scheme() does not error", {
  expect_no_error(mspms:::icelogo_col_scheme())
})

### count_cleavages_per_pos() ###

test_that("count_cleavages_per_pos() generates expected results", {
  test <- tibble::tibble(
    peptide = c("p1", "p2", "p3", "p4"),
    time = c(1, 1, 1, 1),
    condition = c("c1", "c1", "c1", "c1"),
    cleavage_pos = c(1, 1, 1, 1)
  )


  cc <- mspms:::count_cleavages_per_pos(test) %>%
    dplyr::filter(time == 1, condition == "c1", cleavage_pos == 1) %>%
    dplyr::pull(n)


  expect_equal(cc, 4)
})


test_that("count_cleavages_per_pos() generates expected length", {
  seq <- paste0(rep("A", 20),
    collapse = ""
  )

  peptide_library <- tibble::tibble(
    library_id = c("test1"),
    library_match_sequence = seq,
    library_real_sequence = seq
  )

  test <- tibble::tibble(
    peptide = c("p1", "p2", "p3", "p4"),
    time = c(1, 1, 1, 1),
    condition = c("c1", "c1", "c1", "c1"),
    cleavage_pos = c(1, 1, 1, 1)
  )

  cc <- mspms:::count_cleavages_per_pos(test, peptide_library)

  expect_equal(nchar(seq) - 1, max(cc$cleavage_pos))
})
