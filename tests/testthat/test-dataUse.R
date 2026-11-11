test_that(".pdcGQL() refuses to request data without consent", {
    old <- options(CPTACdata.acceptDUA = NULL)
    on.exit(options(old))
    .dataUse[["accepted"]] <- NULL
    # testthat runs non-interactively, so no prompt is shown
    expect_error(.pdcGQL("{ x }"), "Data Use Agreement has not been accepted")
})

test_that("the option records consent without a prompt", {
    old <- options(CPTACdata.acceptDUA = TRUE)
    on.exit(options(old))
    .dataUse[["accepted"]] <- NULL
    expect_true(.checkDataUse())
})
