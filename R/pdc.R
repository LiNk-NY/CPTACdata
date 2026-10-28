# Internal helpers for the Proteomic Data Commons (PDC) GraphQL API.
# Endpoint: https://pdc.cancer.gov/graphql (POST, no auth). Queries that
# return study data need `acceptDUA: true`.

.PDC_URL <- "https://pdc.cancer.gov/graphql"

.pdcGQL <- function(query, timeout = 120) {
    resp <- httr2::request(.PDC_URL) |>
        httr2::req_body_json(list(query = query)) |>
        httr2::req_timeout(timeout) |>
        httr2::req_retry(max_tries = 3L) |>
        httr2::req_error(is_error = function(r) FALSE) |>
        httr2::req_perform()
    body <- httr2::resp_body_json(resp)
    if (httr2::resp_status(resp) >= 400L || !is.null(body[["errors"]])) {
        msgs <- vapply(body[["errors"]],
                       function(e) as.character(e[["message"]]), character(1L))
        stop(sprintf("PDC API error (HTTP %d): %s", httr2::resp_status(resp),
                     paste(msgs, collapse = "; ")), call. = FALSE)
    }
    body[["data"]]
}

.chr <- function(x) if (is.null(x)) NA_character_ else as.character(x)

# Turn a list of GraphQL records into a data.frame, one column per field.
.recordsToDF <- function(records, fields) {
    cols <- lapply(stats::setNames(fields, fields), function(f)
        vapply(records, function(x) .chr(x[[f]]), character(1L)))
    as.data.frame(cols, stringsAsFactors = FALSE)
}

.BIOSPECIMEN_FIELDS <- c(
    "aliquot_submitter_id", "sample_submitter_id", "case_submitter_id",
    "sample_type", "disease_type", "aliquot_status", "pool"
)

# One row per aliquot.
.pdcBiospecimenTable <- local({
    ## only retrieve each study's table once per session, so that a dry run
    ## followed by a real run does not download it twice
    tables <- list()
    function(pdc_study_id) {
        if (is.null(tables[[pdc_study_id]])) {
            q <- sprintf(
                '{ biospecimenPerStudy(pdc_study_id: "%s" acceptDUA: true) { %s } }',
                pdc_study_id, paste(.BIOSPECIMEN_FIELDS, collapse = " ")
            )
            tables[[pdc_study_id]] <<- .recordsToDF(
                .pdcGQL(q)[["biospecimenPerStudy"]], .BIOSPECIMEN_FIELDS
            )
        }
        tables[[pdc_study_id]]
    }
})

# Features x aliquot-label numeric matrix for a study.
.pdcQuantMatrix <- function(pdc_study_id, data_type = "log2_ratio") {
    q <- sprintf(
        '{ quantDataMatrix(pdc_study_id: "%s" data_type: "%s" acceptDUA: true) }',
        pdc_study_id, data_type
    )
    .parseQuantMatrix(.pdcGQL(q, timeout = 600)[["quantDataMatrix"]])
}

# quantDataMatrix returns a list of rows, all character. Row 1 is the header:
# "Gene/Aliquot" (or "Phosphosite/Aliquot") followed by "<uuid>:<label>".
# The uuid is not the PDC aliquot_id; only the label is useful for matching.
# Other rows: feature name followed by values ("NaN" for missing).
.parseQuantMatrix <- function(rows) {
    if (length(rows) < 2L)
        stop("quantDataMatrix returned no data rows.", call. = FALSE)
    header <- vapply(rows[[1L]], .chr, character(1L))
    body <- rows[-1L]
    vals <- suppressWarnings(as.numeric(
        unlist(lapply(body, `[`, -1L), use.names = FALSE)
    ))
    vals[is.nan(vals)] <- NA_real_
    if (length(vals) != length(body) * (length(header) - 1L))
        stop("quantDataMatrix rows have unequal lengths.", call. = FALSE)
    matrix(vals, nrow = length(body), byrow = TRUE, dimnames = list(
        vapply(body, function(r) .chr(r[[1L]]), character(1L)),
        sub("^[^:]+:", "", header[-1L])
    ))
}

#' List PDC studies and their versions
#'
#' Queries \code{studyCatalog} and flattens the result to one row per study
#' version. A PDC study ID (e.g. \code{"PDC000127"}) is stable across versions
#' and resolves to the latest version in other API calls; a version's
#' \code{study_id} (a UUID) pins that version.
#'
#' @param latest_only Logical. If \code{TRUE}, keep only the latest version of
#'   each study.
#' @return A \code{data.frame} with columns \code{pdc_study_id},
#'   \code{study_id}, \code{study_submitter_id}, \code{submitter_id_name},
#'   \code{study_shortname}, \code{study_version}, \code{is_latest_version}.
#' @keywords internal
#' @noRd
.pdcStudyCatalog <- function(latest_only = FALSE) {
    vfields <- c("study_id", "study_submitter_id", "submitter_id_name",
                 "study_shortname", "study_version", "is_latest_version")
    q <- sprintf('{ studyCatalog { pdc_study_id versions { %s } } }',
                 paste(vfields, collapse = " "))
    studies <- .pdcGQL(q, timeout = 60)[["studyCatalog"]]
    ids <- vapply(studies, function(s) .chr(s[["pdc_study_id"]]), character(1L))
    versions <- lapply(studies, `[[`, "versions")
    df <- data.frame(
        pdc_study_id = rep(ids, lengths(versions)),
        .recordsToDF(unlist(versions, recursive = FALSE), vfields),
        stringsAsFactors = FALSE
    )
    # is_latest_version may arrive as a boolean or a "Yes"/"No"-style string
    df$is_latest_version <-
        tolower(df$is_latest_version) %in% c("yes", "true", "1", "y")
    if (latest_only)
        df <- df[df$is_latest_version, , drop = FALSE]
    rownames(df) <- NULL
    df
}
