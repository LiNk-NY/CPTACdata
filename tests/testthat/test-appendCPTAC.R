test_that(".detectTcgaProject() reads admin.disease_code", {
    expect_identical(.detectTcgaProject(mockMAE(brcaPatients, "brca")),
                     "TCGA-BRCA")
    mae <- mockMAE(brcaPatients, "brca")
    MultiAssayExperiment::colData(mae)$admin.disease_code[1L] <- "ov"
    expect_error(.detectTcgaProject(mae), "several projects")
})

test_that(".tcgaSpecimens() drops non-TCGA and pooled specimens", {
    localPDC()
    row <- cptacProjects
    row <- row[row$pdc_study_id == "PDC000173", ]
    spec <- .tcgaSpecimens(row)
    expect_true(all(startsWith(spec$sample_submitter_id, "TCGA-")))
    expect_false(any(tolower(spec$pool) %in% "yes"))
    expect_identical(length(unique(spec$case_submitter_id)), 105L)
})

test_that(".matchQuantColumns() handles both PDC label styles", {
    bio <- data.frame(
        aliquot_submitter_id = c("TCGA-BH-A18U-01A-61-A21W-30",
                                 "TCGA-24-1555-01A-32-2224-29"),
        sample_submitter_id = c("TCGA-BH-A18U-01A", "TCGA-24-1555-01A")
    )
    labels <- c("BH-A18U-01A", "BH-A18U-01A.2", "TCGA-24-1555-01A-32-2224-29",
                "blcdb9-I", "BH-A0E9-01A")
    expect_identical(.matchQuantColumns(labels, bio), c(1L, 1L, 2L, NA, NA))
})

test_that("fraction defaults to a single type, Proteome", {
    localPDC()
    mae <- mockMAE(ovPatients, "ov")
    res <- appendCPTAC(mae, dry_run = TRUE)
    expect_identical(res$experiment,
                     c("CPTAC_JHU_Proteome", "CPTAC_PNNL_Proteome"))
    expect_error(appendCPTAC(mae, fraction = c("Proteome", "Phosphoproteome")),
                 "length 1")
})

test_that("dry_run summarises matches without downloading matrices", {
    localPDC()
    mae <- mockMAE(brcaPatients, "brca")
    res <- appendCPTAC(mae, fraction = "Proteome", dry_run = TRUE)
    expect_s3_class(res, "data.frame")
    expect_identical(res$experiment, "CPTAC_Proteome")
    expect_identical(res$pdc_patients, 105L)
    expect_identical(res$matched_patients, 5L)
    expect_null(gqlCalls[["quant_PDC000173"]])
})

test_that("appendCPTAC() adds a proteome with a valid sampleMap", {
    localPDC()
    mae <- mockMAE(brcaPatients, "brca")
    expect_message(appendCPTAC(mae, fraction = "Proteome", dry_run = TRUE), NA)
    res <- suppressMessages(appendCPTAC(mae, fraction = "Proteome"))
    expect_true(validObject(res))
    expect_identical(names(res), c("RNA", "CPTAC_Proteome"))
    # biospecimen table fetched once across the dry run and real run
    expect_identical(gqlCalls[["bio_PDC000173"]], 1L)

    se <- res[["CPTAC_Proteome"]]
    expect_setequal(colnames(se), paste0(brcaPatients, "-01A"))
    expect_identical(SummarizedExperiment::assayNames(se), "log2_ratio")
    # replicates ".1"/".2" collapse to the first
    expect_identical(se[, "TCGA-C8-A131-01A"]$pdc_column, "C8-A131-01A.1")

    smap <- MultiAssayExperiment::sampleMap(res)
    smap <- smap[smap$assay == "CPTAC_Proteome", ]
    expect_identical(smap$primary, substr(smap$colname, 1L, 12L))

    # re-running does not duplicate the experiment
    expect_message(res2 <- appendCPTAC(res, fraction = "Proteome"),
                   "already present")
    expect_identical(names(res2), names(res))
})

test_that("appendCPTAC() matches full aliquot labels (phosphoproteome)", {
    localPDC()
    mae <- mockMAE(ovPatients, "ov")
    res <- suppressMessages(suppressWarnings(
        appendCPTAC(mae, project = "OV", fraction = "Phosphoproteome")
    ))
    expect_identical(colnames(res[["CPTAC_PNNL_Phosphoproteome"]]),
                     sort(paste0(ovPatients, "-01A")))
})

test_that("appendCPTAC() warns on studies without a matrix", {
    localPDC()
    mae <- mockMAE(ovPatients, "ov")
    expect_warning(
        res <- appendCPTAC(mae, fraction = "Glycoproteome"),
        "no quantification matrix"
    )
    expect_identical(names(res), "RNA")
    expect_identical(ls(gqlCalls), character())
    res <- appendCPTAC(mae, fraction = "Glycoproteome", dry_run = TRUE)
    expect_false(res$quant_available)
    expect_identical(res$matched_patients, NA_integer_)
})

test_that("appendCPTAC() errors for projects without CPTAC data", {
    mae <- mockMAE(brcaPatients, "brca")
    expect_error(appendCPTAC(mae, project = "LUAD"), "No CPTAC Proteome data")
    expect_error(appendCPTAC(matrix(1)), "MultiAssayExperiment")
    expect_error(appendCPTAC(mae, fraction = "Acetylome"), "should be")
})
