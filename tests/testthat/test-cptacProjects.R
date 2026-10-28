test_that("cptacProjects holds the full crosswalk", {
    xw <- cptacProjects
    expect_s3_class(xw, "data.frame")
    expect_identical(nrow(xw), 8L)
    expect_true(all(c("tcga_project", "pdc_study_id", "pdc_name", "experiment",
                      "analytical_fraction", "data_type", "pdc_disease_type",
                      "verified") %in% names(xw)))
    # includes studies without a quantification matrix
    expect_true(anyNA(xw$data_type))
    # experiment names must be unique within a project
    expect_false(anyDuplicated(paste(xw$tcga_project, xw$experiment)) > 0L)
})
