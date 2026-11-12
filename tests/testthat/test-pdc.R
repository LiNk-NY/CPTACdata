test_that(".parseQuantMatrix() reads the header row and converts NaN to NA", {
    rows <- readRDS(test_path("fixtures", "gql_quant_PDC000115.rds"))
    mat <- .parseQuantMatrix(rows[["quantDataMatrix"]])
    expect_true(is.numeric(mat))
    expect_identical(dim(mat), c(5L, 69L))
    expect_identical(rownames(mat)[1L], "NP_000005.2:t118t119")
    # "<uuid>:" prefix removed from column labels
    expect_true(all(startsWith(colnames(mat), "TCGA-")))
    expect_false(any(is.nan(mat)))
    expect_true(is.na(mat[1L, 1L]))
})

test_that(".parseQuantMatrix() keeps annotation columns as rowAnnotation", {
    rows <- list(
        list("Gene/Aliquot", "gene_name", "gene_group",
             "u1:CPT0360520003", "u2:CPT0360520004"),
        list("00e9726d", "VWA3A", "VWA3A", "1.5", "NaN"),
        list("0125e13a", "PGAP2", "PGAP2", "2", "3.25")
    )
    mat <- .parseQuantMatrix(rows)
    expect_identical(colnames(mat), c("CPT0360520003", "CPT0360520004"))
    expect_identical(rownames(mat), c("00e9726d", "0125e13a"))
    expect_identical(mat[, 1L], c(`00e9726d` = 1.5, `0125e13a` = 2))
    expect_true(is.na(mat[1L, 2L]))
    annotation <- attr(mat, "rowAnnotation")
    expect_identical(names(annotation), c("gene_name", "gene_group"))
    expect_identical(annotation$gene_name, c("VWA3A", "PGAP2"))
    # no annotation columns, no attribute
    withoutAnnotation <- lapply(rows, `[`, -(2:3))
    expect_null(attr(.parseQuantMatrix(withoutAnnotation), "rowAnnotation"))
})

test_that(".parseQuantMatrix() errors without data rows", {
    expect_error(.parseQuantMatrix(list(list("Gene/Aliquot"))), "no data rows")
})

test_that(".pdcBiospecimenTable() builds a data.frame and caches it", {
    localPDC()
    bio <- .pdcBiospecimenTable("PDC000173")
    expect_s3_class(bio, "data.frame")
    expect_named(bio, .BIOSPECIMEN_FIELDS)
    expect_identical(nrow(bio), 109L)
    .pdcBiospecimenTable("PDC000173")
    expect_identical(gqlCalls[["bio_PDC000173"]], 1L)
})
