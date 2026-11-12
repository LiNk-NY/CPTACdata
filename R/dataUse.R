# Consent to the PDC Data Use Agreement. Every request to the PDC API goes
# through .pdcGQL(), which calls .checkDataUse() first, so no data are
# requested until the user has agreed.

.PDC_DUA_URL <-
    "https://proteomic.datacommons.cancer.gov/pdc/data-use-guidelines"

.dataUse <- new.env(parent = emptyenv())

.DUA_MESSAGE <- paste0(
    "CPTACdata gets its data from the Proteomic Data Commons (PDC), which\n",
    "is provided under the PDC Data Use Agreement:\n  ", .PDC_DUA_URL, "\n",
    "Read the terms before continuing. Accepting means you agree to them,\n",
    "including any acknowledgement or citation requirements."
)

# How to accept the agreement through the package option, shown both in the
# interactive prompt and in the non-interactive error.
.DUA_OPTION_MESSAGE <- paste0(
    "Accepting at the prompt applies to this R session only. To accept\n",
    "without a prompt, for example in scripts or in every session from your\n",
    ".Rprofile, read the terms and run:\n",
    "  options(CPTACdata.acceptDUA = TRUE)"
)

# Wrapper so that tests can mock interactive().
.isInteractive <- function() {
    interactive()
}

.checkDataUse <- function() {
    if (isTRUE(.dataUse[["accepted"]]) ||
        isTRUE(getOption("CPTACdata.acceptDUA"))) {
        return(invisible(TRUE))
    }
    if (!.isInteractive()) {
        stop(
            "The PDC Data Use Agreement has not been accepted. ",
            .DUA_MESSAGE, "\n",
            .DUA_OPTION_MESSAGE,
            call. = FALSE
        )
    }
    message(.DUA_MESSAGE, "\n\n", .DUA_OPTION_MESSAGE)
    accepted <- BiocBaseUtils::askUserYesNo(
        "Do you accept the PDC Data Use Agreement?"
    )
    if (!accepted) {
        stop(
            "The PDC Data Use Agreement was not accepted; no data were ",
            "requested.",
            call. = FALSE
        )
    }
    .dataUse[["accepted"]] <- TRUE
    invisible(TRUE)
}
