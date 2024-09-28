test_that("prepare_pd() ", {
    peptide_groups_filepath <- system.file("extdata/proteome_discoverer_PeptideGroups.txt",
        package = "mspms"
    )
    colData_filepath <- system.file("extdata/proteome_discover_colData.csv",
        package = "mspms"
    )
    expect_s4_class(
        mspms::prepare_pd(
            peptide_groups_filepath,
            colData_filepath
        ),
        "QFeatures"
    )
})
