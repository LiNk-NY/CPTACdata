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

test_that("the prompt tells the user about the CPTACdata.acceptDUA option", {
    old <- options(CPTACdata.acceptDUA = NULL)
    on.exit(options(old))
    .dataUse[["accepted"]] <- NULL
    on.exit(.dataUse[["accepted"]] <- NULL, add = TRUE)
    local_mocked_bindings(.isInteractive = function() TRUE)
    local_mocked_bindings(
        askUserYesNo = function(...) TRUE,
        .package = "BiocBaseUtils"
    )
    expect_message(
        expect_true(.checkDataUse()),
        "options(CPTACdata.acceptDUA = TRUE)",
        fixed = TRUE
    )
    expect_true(.dataUse[["accepted"]])
})

test_that("declining at the prompt stops without recording consent", {
    old <- options(CPTACdata.acceptDUA = NULL)
    on.exit(options(old))
    .dataUse[["accepted"]] <- NULL
    local_mocked_bindings(.isInteractive = function() TRUE)
    local_mocked_bindings(
        askUserYesNo = function(...) FALSE,
        .package = "BiocBaseUtils"
    )
    expect_error(
        suppressMessages(.checkDataUse()),
        "was not accepted"
    )
    expect_null(.dataUse[["accepted"]])
})

test_that("the non-interactive error names the option", {
    old <- options(CPTACdata.acceptDUA = NULL)
    on.exit(options(old))
    .dataUse[["accepted"]] <- NULL
    expect_error(
        .checkDataUse(),
        "options(CPTACdata.acceptDUA = TRUE)",
        fixed = TRUE
    )
})
