## Build the `cptacCohorts` table of CPTAC cohorts and save it as exported
## package data (data/cptacCohorts.rda).
##
## Run from the package root: Rscript data-raw/cptacCohorts.R
## Needs internet access to pdc.cancer.gov, and accepts the PDC Data Use
## Agreement for the queries it makes.
##
## Studies come from the PDC projects "CPTAC3 Discovery and Confirmatory"
## and "CPTAC-2" (the Prospective collections). CompRef studies (reference
## samples only) and the studies in `excluded` are left out. In CPTAC-3
## study names are "<cohort> - <experiment>"; the CPTAC-2 names have no
## separator, so their cohort and experiment are curated in `prospective`.
## The script stops if a study name does not fit these rules or if the
## finished table fails the checks in `checkCohorts()`.
options(CPTACdata.acceptDUA = TRUE)
library(tibble)
library(dplyr)
devtools::load_all(quiet = TRUE)

.COHORT_PROJECTS <- c("CPTAC3 Discovery and Confirmatory", "CPTAC-2")

## Curated: studies in those projects that are not patient cohorts.
excluded <- c(
    ## method comparison on PDAC tissue, not a cohort
    "PDC000393"
)

## Curated: cohort and experiment names of the CPTAC-2 Prospective studies.
prospective <- tribble(
    ~pdc_study_id , ~cohort                , ~experiment                ,
    "PDC000120"   , "Prospective Breast"   , "BI_Proteome"              ,
    "PDC000121"   , "Prospective Breast"   , "BI_Phosphoproteome"       ,
    "PDC000239"   , "Prospective Breast"   , "BI_Acetylome"             ,
    "PDC000109"   , "Prospective Colon"    , "VU_Proteome"              ,
    "PDC000116"   , "Prospective Colon"    , "PNNL_Proteome"            ,
    "PDC000117"   , "Prospective Colon"    , "PNNL_Phosphoproteome"     ,
    "PDC000110"   , "Prospective Ovarian"  , "JHU_Proteome"             ,
    "PDC000250"   , "Prospective Ovarian"  , "JHU_Intact_Glycoproteome" ,
    "PDC000251"   , "Prospective Ovarian"  , "JHU_Glycosite_Peptides"   ,
    "PDC000118"   , "Prospective Ovarian"  , "PNNL_Proteome"            ,
    "PDC000119"   , "Prospective Ovarian"  , "PNNL_Phosphoproteome"
)

## quantDataMatrix data types tried, in order; the first that the PDC has
## is recorded.
.DATA_TYPES <- c(
    "log2_ratio",
    "unshared_log2_ratio",
    "precursor_area",
    "unshared_precursor_area",
    "spectral_count",
    "unshared_spectral_count"
)

## All studies of the cohort projects, from allPrograms.
pdcProjectStudies <- function(projects) {
    progs <- .pdcGQL(paste0(
        '{ allPrograms(acceptDUA: true) { projects { project_submitter_id ',
        'studies { pdc_study_id submitter_id_name analytical_fraction ',
        'experiment_type } } } }'
    ))
    rows <- list()
    for (program in progs[["allPrograms"]]) {
        for (project in program[["projects"]]) {
            if (!project[["project_submitter_id"]] %in% projects) {
                next
            }
            for (study in project[["studies"]]) {
                rows[[length(rows) + 1L]] <- tibble(
                    pdc_project = project[["project_submitter_id"]],
                    pdc_study_id = study[["pdc_study_id"]],
                    pdc_name = study[["submitter_id_name"]],
                    analytical_fraction = .toCharacter(
                        study[["analytical_fraction"]]
                    ),
                    experiment_type = .toCharacter(study[["experiment_type"]])
                )
            }
        }
    }
    distinct(bind_rows(rows))
}

## Column labels of the study's quantDataMatrix of this type, or NULL if the
## PDC has none. Only the header row is read, so whole matrices are not
## downloaded.
matrixLabels <- function(pdc_study_id, data_type) {
    query <- sprintf(
        paste0(
            '{ quantDataMatrix(pdc_study_id: "%s" data_type: "%s" ',
            'acceptDUA: true) }'
        ),
        pdc_study_id,
        data_type
    )
    connection <- httr2::request(.PDC_URL) |>
        httr2::req_body_json(list(query = query)) |>
        httr2::req_timeout(600) |>
        httr2::req_retry(
            max_tries = 3L,
            is_transient = function(response) {
                httr2::resp_status(response) %in% .PDC_TRANSIENT_STATUS
            }
        ) |>
        httr2::req_error(is_error = function(response) FALSE) |>
        httr2::req_perform_connection(blocking = TRUE)
    on.exit(close(connection))
    start <- ""
    repeat {
        chunk <- httr2::resp_stream_raw(connection, kb = 16)
        start <- paste0(start, rawToChar(chunk))
        if (!length(chunk) || grepl('\\[\\[[^]]*\\]', start)) {
            break
        }
    }
    if (!grepl('"quantDataMatrix":\\[\\[', start)) {
        return(NULL)
    }
    header <- regmatches(start, regexpr('\\[\\[[^]]*\\]', start))
    labels <- regmatches(header, gregexpr('"[^"]*"', header))[[1L]]
    labels <- gsub('"', "", labels)[-1L]
    sub("^[^:]+:", "", labels)
}

## Matrix labels of reference channels whose label does not match their
## biospecimen record, e.g. "POOL.1", "WU-POOL", "none", or CompRef pools
## such as "LUAD_Global_CR_pool_1.1" and "JHU_HNSCC_CR.2" (recorded as
## "LUAD Global CR (pool 1)" and "JHU HNSCC CR").
.REFERENCE_LABELS <- "(?i)pool|^none(\\.[0-9]+)?$|_CR([._]|$)"

## The first data type in .DATA_TYPES with a matrix, and the share of its
## sample labels that match an aliquot of the study. Annotation columns
## (.QUANT_ANNOTATION_COLUMNS) and reference channels are not sample labels.
matrixDataType <- function(pdc_study_id) {
    for (data_type in .DATA_TYPES) {
        labels <- matrixLabels(pdc_study_id, data_type)
        if (!is.null(labels)) {
            labels <- labels[
                !labels %in% .QUANT_ANNOTATION_COLUMNS &
                    !grepl(.REFERENCE_LABELS, labels, perl = TRUE)
            ]
            biospecimens <- .pdcBiospecimenTable(pdc_study_id)
            matched <- !is.na(.matchAliquotLabels(labels, biospecimens))
            return(tibble(data_type = data_type, labels_matched = mean(matched)))
        }
    }
    tibble(data_type = NA_character_, labels_matched = NA_real_)
}

## "CPTAC CCRCC Discovery Study - DIA Proteome" -> cohort
## "CPTAC CCRCC Discovery Study", experiment "DIA_Proteome".
splitStudyName <- function(pdc_name) {
    parts <- strsplit(pdc_name, " - ", fixed = TRUE)
    tibble(
        cohort = vapply(parts, `[`, character(1L), 1L),
        experiment = vapply(parts, function(part) {
            name <- paste(part[-1L], collapse = " ")
            name <- gsub("[^A-Za-z0-9]+", "_", name)
            gsub("^_|_$", "", name)
        }, character(1L))
    )
}

## ------------------------------------------------------------- checks --

## Invariants of the finished table. Returns the failed checks as a
## character vector (empty if all pass); the build stops if there are any.
checkCohorts <- function(cohorts) {
    failures <- character()
    check <- function(ok, message) {
        if (!isTRUE(ok)) {
            failures <<- c(failures, message)
        }
    }
    check(
        !anyNA(cohorts[c("cohort", "pdc_study_id", "experiment")]) &&
            all(nzchar(cohorts[["experiment"]])),
        "missing cohort, pdc_study_id or experiment"
    )
    check(
        all(grepl("^PDC[0-9]{6}$", cohorts[["pdc_study_id"]])) &&
            !anyDuplicated(cohorts[["pdc_study_id"]]),
        "malformed or repeated pdc_study_id"
    )
    check(
        !anyDuplicated(paste(cohorts[["cohort"]], cohorts[["experiment"]])),
        "experiment names repeat within a cohort"
    )
    check(
        all(grepl("^[A-Za-z][A-Za-z0-9_]*$", cohorts[["experiment"]])),
        "experiment names are not syntactic"
    )
    check(
        !any(grepl("CompRef", cohorts[["pdc_name"]])),
        "CompRef studies included"
    )
    check(
        all(cohorts[["data_type"]] %in% c(.DATA_TYPES, NA)),
        "unexpected data_type"
    )
    check(
        all(tapply(!is.na(cohorts[["data_type"]]), cohorts[["cohort"]], any)),
        "a cohort has no study with a quantification matrix"
    )
    check(
        all(cohorts[["labels_matched"]] == 1, na.rm = TRUE),
        "matrix labels that are not aliquot_submitter_ids"
    )
    failures
}

## --------------------------------------------------------------- build --

studies <- pdcProjectStudies(.COHORT_PROJECTS) |>
    filter(
        !grepl("CompRef", pdc_name),
        !pdc_study_id %in% excluded
    )

isProspective <- studies[["pdc_project"]] == "CPTAC-2"
unknownProspective <- setdiff(
    studies[["pdc_study_id"]][isProspective],
    prospective[["pdc_study_id"]]
)
if (length(unknownProspective)) {
    stop(
        "CPTAC-2 studies missing from `prospective`: ",
        paste(unknownProspective, collapse = ", ")
    )
}
cptac3 <- studies[!isProspective, ]
unsplit <- cptac3[["pdc_study_id"]][!grepl(" - ", cptac3[["pdc_name"]])]
if (length(unsplit)) {
    stop(
        "CPTAC-3 study names without ' - ': ",
        paste(unsplit, collapse = ", ")
    )
}
names3 <- bind_cols(
    cptac3["pdc_study_id"],
    splitStudyName(cptac3[["pdc_name"]])
)

cptacCohorts <- studies |>
    inner_join(bind_rows(prospective, names3), by = "pdc_study_id")
cptacCohorts <- bind_cols(
    cptacCohorts,
    bind_rows(lapply(cptacCohorts[["pdc_study_id"]], matrixDataType))
)

## Cohorts without any quantification matrix cannot be downloaded
noMatrix <- cptacCohorts |>
    group_by(cohort) |>
    filter(all(is.na(data_type))) |>
    ungroup()
if (nrow(noMatrix)) {
    message(
        "Dropping cohorts without a quantification matrix: ",
        paste(unique(noMatrix[["cohort"]]), collapse = ", ")
    )
}
cptacCohorts <- cptacCohorts |>
    filter(!cohort %in% noMatrix[["cohort"]])

failures <- checkCohorts(cptacCohorts)
if (length(failures)) {
    print(as.data.frame(cptacCohorts[
        cptacCohorts[["labels_matched"]] %in% setdiff(
            cptacCohorts[["labels_matched"]], c(1, NA)
        ),
        c("pdc_study_id", "pdc_name", "data_type", "labels_matched")
    ]))
    stop("cptacCohorts failed checks:\n  ", paste(failures, collapse = "\n  "))
}

cptacCohorts <- cptacCohorts |>
    select(
        cohort, experiment, pdc_study_id, pdc_name, pdc_project,
        analytical_fraction, experiment_type, data_type
    ) |>
    arrange(cohort, experiment) |>
    as.data.frame()
save(cptacCohorts, file = "data/cptacCohorts.rda", compress = "xz")
