test_that("mspms_tidy() returns a tibble with the expected colnames", {
    expected_colData_names <- c("quantCols", "time", "condition", "group")
    expected_rowData_names <- c(
        "peptide",
        "peptide_type",
        "library_match_sequence",
        "library_real_sequence",
        "library_id",
        "cleavage_seq",
        "cleavage_pos"
    )
    expected_value_name <- "peptides_norm"
    expected_names <- c(
        expected_colData_names, expected_rowData_names,
        expected_value_name
    )
    tidy <- mspms_tidy(mspms::processed_qf)
    colnames <- colnames(tidy)
    n_missing <- sum(expected_names %!in% colnames)
    expect_equal(n_missing, 0)
})
