ccrcc <- "CPTAC CCRCC Discovery Study"

test_that("cptacCohorts lists the CCRCC discovery studies", {
    cohorts <- .loadData("cptacCohorts")
    expect_s3_class(cohorts, "data.frame")
    expect_false(anyDuplicated(cohorts$pdc_study_id) > 0L)
    expect_false(any(grepl("CompRef", cohorts$pdc_name)))
    rows <- cohorts[cohorts$cohort == ccrcc, ]
    expect_identical(
        rows$pdc_study_id[match(c("Proteome", "Phosphoproteome"),
                                rows$experiment)],
        c("PDC000127", "PDC000128")
    )
})

test_that(".matchAliquotLabels() matches labels with replicate suffixes", {
    bio <- data.frame(aliquot_submitter_id = c("CPT0026410003",
                                               "CPT0002370001"))
    labels <- c("CPT0026410003", "CPT0002370001.2", "QC5", "NCI7-3")
    expect_identical(.matchAliquotLabels(labels, bio), c(1L, 2L, NA, NA))
})

test_that(".cohortSpecimens() drops reference channels", {
    localPDC()
    row <- data.frame(pdc_study_id = "PDC000127")
    spec <- .cohortSpecimens(row)
    expect_true(all(grepl("^C3[LN]-", spec$case_submitter_id)))
    expect_false(any(spec$disease_type %in% "Other"))
})

test_that("getCPTAC() dry run summarises the cohort", {
    localPDC()
    res <- getCPTAC(ccrcc, c("Proteome", "Phosphoproteome"), dry_run = TRUE)
    expect_identical(res$experiment, c("Proteome", "Phosphoproteome"))
    expect_true(all(res$quant_available))
    expect_true(all(res$patients > 100L))
    expect_null(gqlCalls[["quant_PDC000127"]])
})

test_that("getCPTAC() combines studies in a MultiAssayExperiment", {
    localPDC()
    res <- suppressMessages(
        getCPTAC(ccrcc, c("Proteome", "Phosphoproteome"))
    )
    expect_s4_class(res, "MultiAssayExperiment")
    expect_true(validObject(res))
    expect_identical(names(res), c("Proteome", "Phosphoproteome"))

    se <- res[["Proteome"]]
    expect_identical(SummarizedExperiment::assayNames(se), "log2_ratio")
    # columns keep the PDC labels, which are aliquot_submitter_ids
    expect_identical(colnames(se), se$pdc_column)
    expect_true(all(startsWith(colnames(se), "CPT")))
    expect_identical(S4Vectors::metadata(se)[["pdc_study_id"]], "PDC000127")

    smap <- MultiAssayExperiment::sampleMap(res)
    expect_true(all(grepl("^C3[LN]-", smap$primary)))
    cd <- MultiAssayExperiment::colData(res)
    expect_identical(rownames(cd), cd$case_submitter_id)
    expect_true(all(c("gender", "primary_diagnosis", "tumor_grade") %in%
                    names(cd)))
    expect_false(anyNA(cd$disease_type))
    # tumor and normal samples of a patient map to the same primary
    expect_true(any(table(smap$primary[smap$assay == "Proteome"]) > 1L))
})

test_that("getCPTAC() errors on unknown cohorts and experiments", {
    expect_error(getCPTAC("Not a cohort"), "not a cohort")
    expect_error(getCPTAC(ccrcc, "Lipidome_X"), "Not experiments")
    expect_error(getCPTAC(ccrcc, replicates = "mean"), "should be")
})
