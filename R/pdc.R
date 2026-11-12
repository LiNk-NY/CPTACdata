# Internal helpers for the Proteomic Data Commons (PDC) GraphQL API.
# Endpoint: https://pdc.cancer.gov/graphql (POST, no auth). Queries that
# return study data need `acceptDUA: true`.

.PDC_URL <- "https://pdc.cancer.gov/graphql"

# HTTP statuses retried as temporary server or gateway failures.
.PDC_TRANSIENT_STATUS <- c(429L, 500L, 502L, 503L, 504L)

.pdcGQL <- function(query, timeout = 120) {
    .checkDataUse()
    response <- httr2::request(.PDC_URL) |>
        httr2::req_body_json(list(query = query)) |>
        httr2::req_timeout(timeout) |>
        httr2::req_retry(
            max_tries = 3L,
            is_transient = function(response) {
                httr2::resp_status(response) %in% .PDC_TRANSIENT_STATUS
            }
        ) |>
        httr2::req_error(is_error = function(response) FALSE) |>
        httr2::req_perform()
    if (!identical(httr2::resp_content_type(response), "application/json")) {
        stop(
            "PDC API request failed (HTTP ",
            httr2::resp_status(response), "): the response was not JSON.",
            call. = FALSE
        )
    }
    body <- httr2::resp_body_json(response)
    if (httr2::resp_status(response) >= 400L || !is.null(body[["errors"]])) {
        messages <- vapply(
            body[["errors"]],
            function(error) as.character(error[["message"]]),
            character(1L)
        )
        stop(
            "PDC API request failed (HTTP ",
            httr2::resp_status(response), "): ",
            paste(messages, collapse = "; "),
            call. = FALSE
        )
    }
    body[["data"]]
}

.toCharacter <- function(value) {
    if (is.null(value)) NA_character_ else as.character(value)
}

# Turn a list of GraphQL records into a data.frame, one column per field.
.recordsToDF <- function(records, fields) {
    columns <- lapply(stats::setNames(fields, fields), function(field) {
        vapply(
            records,
            function(record) .toCharacter(record[[field]]),
            character(1L)
        )
    })
    as.data.frame(columns, stringsAsFactors = FALSE)
}

.BIOSPECIMEN_FIELDS <- c(
    "aliquot_submitter_id",
    "sample_submitter_id",
    "case_submitter_id",
    "sample_type",
    "disease_type",
    "aliquot_status",
    "pool"
)

# Session cache for values fetched from the PDC: biospecimen tables and the
# study catalog. Each is fetched once, so a dry run followed by a real run
# does not download it twice.
.pdcCache <- new.env(parent = emptyenv())

# One row per aliquot.
.pdcBiospecimenTable <- function(pdc_study_id) {
    key <- paste0("biospecimen_", pdc_study_id)
    if (is.null(.pdcCache[[key]])) {
        query <- sprintf(
            paste0(
                '{ biospecimenPerStudy(pdc_study_id: "%s" ',
                'acceptDUA: true) { %s } }'
            ),
            pdc_study_id,
            paste(.BIOSPECIMEN_FIELDS, collapse = " ")
        )
        .pdcCache[[key]] <- .recordsToDF(
            .pdcGQL(query)[["biospecimenPerStudy"]],
            .BIOSPECIMEN_FIELDS
        )
    }
    .pdcCache[[key]]
}

# Features x aliquot-label numeric matrix for a study.
.pdcQuantMatrix <- function(pdc_study_id, data_type = "log2_ratio") {
    query <- sprintf(
        paste0(
            '{ quantDataMatrix(pdc_study_id: "%s" data_type: "%s" ',
            'acceptDUA: true) }'
        ),
        pdc_study_id,
        data_type
    )
    .parseQuantMatrix(.pdcGQL(query, timeout = 600)[["quantDataMatrix"]])
}

# Matrix columns that annotate the features rather than hold sample values.
# Newer studies identify rows by a PDC gene ID and give the gene symbol in
# these columns.
.QUANT_ANNOTATION_COLUMNS <- c("gene_name", "gene_group")

# quantDataMatrix returns a list of rows, all character. Row 1 is the header:
# "Gene/Aliquot" (or "Phosphosite/Aliquot"), in newer studies the annotation
# columns in .QUANT_ANNOTATION_COLUMNS, then "<uuid>:<label>". The uuid is
# not the PDC aliquot_id; only the label is useful for matching. Other rows:
# feature ID, annotation values, then sample values ("NaN" for missing).
# Returns the numeric matrix, with the annotation columns, if any, as a
# data.frame in the "rowAnnotation" attribute.
.parseQuantMatrix <- function(rows) {
    if (length(rows) < 2L) {
        stop("quantDataMatrix returned no data rows.", call. = FALSE)
    }
    header <- vapply(rows[[1L]], .toCharacter, character(1L))
    labels <- sub("^[^:]+:", "", header[-1L])
    body <- rows[-1L]
    cells <- unlist(lapply(body, `[`, -1L), use.names = FALSE)
    if (length(cells) != length(body) * length(labels)) {
        stop("quantDataMatrix rows have unequal lengths.", call. = FALSE)
    }
    cells <- matrix(cells, nrow = length(body), byrow = TRUE)
    featureIds <- vapply(
        body,
        function(row) .toCharacter(row[[1L]]),
        character(1L)
    )
    isAnnotation <- labels %in% .QUANT_ANNOTATION_COLUMNS
    # "NaN" converts to NaN without a warning, so only real non-numeric text
    # would warn
    values <- as.numeric(cells[, !isAnnotation, drop = FALSE])
    values[is.nan(values)] <- NA_real_
    quantMatrix <- matrix(
        values,
        nrow = length(body),
        dimnames = list(featureIds, labels[!isAnnotation])
    )
    if (any(isAnnotation)) {
        attr(quantMatrix, "rowAnnotation") <- as.data.frame(
            cells[, isAnnotation, drop = FALSE],
            row.names = featureIds,
            stringsAsFactors = FALSE
        ) |>
            stats::setNames(labels[isAnnotation])
    }
    quantMatrix
}

#' List PDC studies and their versions
#'
#' Queries `studyCatalog` and flattens the result to one row per study version.
#' A PDC study ID (e.g. `"PDC000127"`) is stable across versions and resolves
#' to the latest version in other API calls; a version's `study_id` (a UUID)
#' pins that version.
#'
#' @param latest_only Logical. If `TRUE`, keep only the latest version of each
#'   study.
#' @return A `data.frame` with columns `pdc_study_id`, `study_id`,
#'   `study_submitter_id`, `submitter_id_name`, `study_shortname`,
#'   `study_version`, `is_latest_version`.
#' @keywords internal
#' @noRd
.pdcStudyCatalog <- function(latest_only = FALSE) {
    versionFields <- c(
        "study_id",
        "study_submitter_id",
        "submitter_id_name",
        "study_shortname",
        "study_version",
        "is_latest_version"
    )
    query <- sprintf(
        '{ studyCatalog { pdc_study_id versions { %s } } }',
        paste(versionFields, collapse = " ")
    )
    studies <- .pdcGQL(query, timeout = 60)[["studyCatalog"]]
    studyIds <- vapply(
        studies,
        function(study) .toCharacter(study[["pdc_study_id"]]),
        character(1L)
    )
    versions <- lapply(studies, `[[`, "versions")
    catalog <- data.frame(
        pdc_study_id = rep(studyIds, lengths(versions)),
        .recordsToDF(unlist(versions, recursive = FALSE), versionFields),
        stringsAsFactors = FALSE
    )
    # is_latest_version may arrive as a boolean or a "Yes"/"No"-style string
    catalog[["is_latest_version"]] <-
        tolower(catalog[["is_latest_version"]]) %in% c("yes", "true", "1", "y")
    if (latest_only) {
        catalog <- catalog[catalog[["is_latest_version"]], , drop = FALSE]
    }
    rownames(catalog) <- NULL
    catalog
}

# Version of a study that the PDC resolves its pdc_study_id to.
.pdcStudyVersion <- function(pdc_study_id) {
    if (is.null(.pdcCache[["studyCatalog"]])) {
        .pdcCache[["studyCatalog"]] <- .pdcStudyCatalog(latest_only = TRUE)
    }
    catalog <- .pdcCache[["studyCatalog"]]
    match <- catalog[catalog[["pdc_study_id"]] == pdc_study_id, ,
        drop = FALSE]
    if (nrow(match) != 1L) {
        stop(
            "PDC catalog has ", nrow(match), " latest versions of '",
            pdc_study_id, "'.",
            call. = FALSE
        )
    }
    list(
        pdc_study_version = match[["study_version"]],
        pdc_study_uuid = match[["study_id"]]
    )
}

# ------------------------------------------------------------- clinical --

.DEMOGRAPHIC_FIELDS <- c(
    "gender",
    "race",
    "ethnicity",
    "vital_status",
    "days_to_birth",
    "days_to_death",
    "cause_of_death"
)

.DIAGNOSIS_FIELDS <- c(
    "primary_diagnosis",
    "morphology",
    "tumor_grade",
    "tumor_stage",
    "ajcc_pathologic_stage",
    "age_at_diagnosis",
    "tissue_or_organ_of_origin",
    "site_of_resection_or_biopsy",
    "days_to_last_follow_up",
    "last_known_disease_status"
)

# All records of a paginated per-study case query, e.g.
# paginatedCaseDemographicsPerStudy, fetched `limit` cases at a time.
.pdcPaginatedCases <- function(
    pdc_study_id, queryName, listName, fields, limit = 500L
) {
    records <- list()
    offset <- 0L
    repeat {
        query <- sprintf(
            paste0(
                '{ %s(pdc_study_id: "%s" offset: %d limit: %d ',
                'acceptDUA: true) { total %s { %s } } }'
            ),
            queryName, pdc_study_id, offset, limit, listName, fields
        )
        page <- .pdcGQL(query)[[queryName]]
        records <- c(records, page[[listName]])
        offset <- offset + limit
        if (offset >= page[["total"]] || !length(page[[listName]])) {
            break
        }
    }
    records
}

# One row per case from the nested demographics or diagnoses records. A case
# with several records gets their distinct values joined with "; ".
.flattenCaseRecords <- function(records, nestedName, fields) {
    rows <- lapply(records, function(record) {
        nested <- .recordsToDF(record[[nestedName]], fields)
        values <- vapply(fields, function(field) {
            distinct <- unique(stats::na.omit(nested[[field]]))
            if (length(distinct)) {
                paste(distinct, collapse = "; ")
            } else {
                NA_character_
            }
        }, character(1L))
        data.frame(
            case_submitter_id = .toCharacter(record[["case_submitter_id"]]),
            as.list(values),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    })
    if (!length(rows)) {
        empty <- rep(list(character()), length(fields) + 1L)
        names(empty) <- c("case_submitter_id", fields)
        return(as.data.frame(empty, stringsAsFactors = FALSE))
    }
    do.call(rbind, rows)
}

# Clinical data for the cases of a study: one row per case, with
# case_submitter_id, disease_type, demographics and diagnoses.
.pdcClinical <- function(pdc_study_id) {
    key <- paste0("clinical_", pdc_study_id)
    if (is.null(.pdcCache[[key]])) {
        demographics <- .pdcPaginatedCases(
            pdc_study_id,
            "paginatedCaseDemographicsPerStudy",
            "caseDemographicsPerStudy",
            sprintf(
                "case_submitter_id disease_type demographics { %s }",
                paste(.DEMOGRAPHIC_FIELDS, collapse = " ")
            )
        )
        diagnoses <- .pdcPaginatedCases(
            pdc_study_id,
            "paginatedCaseDiagnosesPerStudy",
            "caseDiagnosesPerStudy",
            sprintf(
                "case_submitter_id diagnoses { %s }",
                paste(.DIAGNOSIS_FIELDS, collapse = " ")
            )
        )
        cases <- .recordsToDF(
            demographics,
            c("case_submitter_id", "disease_type")
        )
        clinical <- merge(
            cases,
            .flattenCaseRecords(
                demographics, "demographics", .DEMOGRAPHIC_FIELDS
            ),
            by = "case_submitter_id",
            all.x = TRUE,
            sort = FALSE
        )
        .pdcCache[[key]] <- merge(
            clinical,
            .flattenCaseRecords(diagnoses, "diagnoses", .DIAGNOSIS_FIELDS),
            by = "case_submitter_id",
            all.x = TRUE,
            sort = FALSE
        )
    }
    .pdcCache[[key]]
}
