test_that("prepare_fragpipe() results in a QFeatures object", {
    prepared_data <- prepare_fragpipe(
        combined_peptide_filepath <- system.file(
            "extdata/fragpipe_combined_peptide.tsv",
            package = "mspms"
        ),
        colData_filepath <- system.file("extdata/colData.csv", package = "mspms")
    )

    expect_s4_class(prepared_data, "QFeatures")
})

test_that("prepare_fragpipe() produces error when
          peaks data is accidently loaded", {
    expect_error(prepare_fragpipe(
        combined_peptide_filepath <- system.file(
            "extdata/peaks_protein-peptides-lfq.csv",
            package = "mspms"
        ),
        colData_filepath <- system.file("extdata/colData.csv", package = "mspms")
    ))
})
