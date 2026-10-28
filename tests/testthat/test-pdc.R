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
